#include <cstdlib>
#include <cmath>

#include "common.hh"
#include "geometry.hh"
#include "region.hh"

/** Initializes the geometry class, setting up the processor grid to minimize
 * the surface area between processors.
 * \param[in] (m_,n_,o_) the global grid dimensions.
 * \param[in] (x_prd_,y_prd_,z_prd_) the global periodicities. */
geometry::geometry(int m_,int n_,int o_,bool x_prd_,bool y_prd_,bool z_prd_)
	: m(m_), n(n_), o(o_), x_prd(x_prd_), y_prd(y_prd_), z_prd(z_prd_),
	ranks(NULL) {

	// Determine the total number of processors
	MPI_Comm_size(MPI_COMM_WORLD,&procs);

	// Find the optimal grid of processors to minimize the total communication
	long msize=long(m)*long(n)*long(o);
	search_optimal(msize);
	if(msize==m*n*o) p_fatal_error("No valid processor decomposition",PGMG_MATH_ERROR);

	// Create the Cartesian communicator, and set the grid extent and rank.
	create_cart_communicator(MPI_COMM_WORLD);
}

/** Initializes the geometry class from a region class that represents part of
 * a grid in a multigrid hierarchy. It uses the same global grid dimensions and
 * topology but thins out the number of processors used by a given factor.
 * \param[in] p a pointer to the region class.
 * \param[in] pprocs the number of processors involved in the region class.
 * \param[in] thin the thinning factor. */
geometry::geometry(region* p,int pprocs,int thin) : m(p->m), n(p->n), o(p->o), x_prd(p->x_prd),
	y_prd(p->y_prd), z_prd(p->z_prd) {

	// Reduce the number of processors and determine the optimal processor
	// decomposition
	procs=pprocs/thin;
	if(procs<1) procs=1;
	long msize=long(m)*long(n)*long(o);
	search_optimal(msize);
	if(msize==m*n*o) p_fatal_error("No valid processor decomposition",PGMG_MATH_ERROR);

	// Create group of processors by thinning out the previous set
	MPI_Comm new_comm;
	MPI_Comm_split(p->cart,p->rank<procs?0:MPI_UNDEFINED,p->rank,&new_comm);

	// Make new communicator. If this processor is in the communicator,
	// then set up all of the constants.
	if(new_comm==MPI_COMM_NULL) rank=-1;
	else {
		create_cart_communicator(new_comm);
		MPI_Comm_free(&new_comm);
	}

	// Set up the global table of processor ranks. This is needed for the
	// parent regions to find the right children to talk to, since some
	// parents are not in the children's communicator.
	ranks=new int[mp*np*op];
	if(rank==0) {
		int *rp=ranks,q[3],&i=*q,&j=q[1],&k=q[2];
		for(k=0;k<op;k++) for(j=0;j<np;j++) for(i=0;i<mp;i++,rp++)
			MPI_Cart_rank(cart,q,rp);
	}
	MPI_Bcast(ranks,mp*np*op,MPI_INT,0,p->cart);
}

/**
 * Given a geometry with dimensions (m,n,o), initializes a new geometry
 * with the same processor decomposition. In cartesian directions that
 * have non-periodic boundaries, the new geometry will have an additional
 * grid point on the last plane of processors.
 */
// TODO don't overload geometry(const geometry &gm)?
geometry::geometry(const geometry &gm) :
	m(gm.x_prd?gm.m:gm.m+1),
	n(gm.y_prd?gm.n:gm.n+1),
	o(gm.z_prd?gm.o:gm.o+1),
	x_prd(gm.x_prd),y_prd(gm.y_prd),z_prd(gm.z_prd),
	rank(gm.rank),procs(gm.procs),
	mp(gm.mp),np(gm.np),op(gm.op),
	ip(gm.ip),jp(gm.jp),kp(gm.kp),
	ai(gm.ai),aj(gm.aj),ak(gm.ak),
	bi(gm.x_prd||gm.ip<(gm.mp-1)?gm.bi:gm.bi+1),
	bj(gm.y_prd||gm.jp<(gm.np-1)?gm.bj:gm.bj+1),
	bk(gm.z_prd||gm.kp<(gm.op-1)?gm.bk:gm.bk+1),
	sm(bi-ai),sn(bj-aj),so(bk-ak),
	smn(sm*sn),smno(smn*so),ranks(NULL) {

	// the above copies all data over to the new geometry, incrementing
	// b(i/j/k) and s(m/n/o) if (x/y/z)_prd is false and
	// (i/j/k)p = (m/n/o)p-1. To finish, duplicate the communicator
	MPI_Comm_dup(gm.cart,&cart);
}

/** The geometry destructor frees the dynamically allocated global table of
 * ranks, if present. */
geometry::~geometry() {
	if(ranks!=NULL) delete [] ranks;
}

/** Frees the memory used in the Cartesian communicator. */
void geometry::free() {
	if(rank!=-1) MPI_Comm_free(&cart);
}

/** Searches for an optimal processor decomposition, to minimize the surface area between
 * processors.
 * \param[in,out] msize the current best surface area, which is replaced by a
 *			lower surface area if such a decomposition is found. */
void geometry::search_optimal(long &msize) {
	int i,j,k,p2;
	long tsize;
	for(k=1;k<=procs;k++) {
		if(procs%k!=0) continue;
		p2=procs/k;
		for(j=1;j<=p2;j++) {
			if(p2%j!=0) continue;
			i=p2/j;

			// Calculate the surface area betweeen the processors
			// for this particular processor decomposition
			tsize=m*n*(z_prd?k:k-1)+m*o*(y_prd?j:j-1)+n*o*(x_prd?i:i-1);

			// If this surface area is smaller than the previous
			// smallest known, then remember this processor
			// decomposition
			if(tsize<msize) {
				msize=tsize;
				mp=i;np=j;op=k;
			}
		}
	}
}

/** Sets up and stores the Cartesian communicator and associated constants.
 * \param[in] mcom the MPI communicator to use. */
void geometry::create_cart_communicator(MPI_Comm mcom) {

	// Create communicator
	int dims[3],pers[3],coor[3];
	*dims=mp;dims[1]=np;dims[2]=op;
	*pers=x_prd;pers[1]=y_prd;pers[2]=z_prd;
	MPI_Cart_create(mcom,3,dims,pers,1,&cart);

	// Set this processor's rank and position within the grid
	MPI_Comm_rank(cart,&rank);
	MPI_Cart_coords(cart,rank,3,coor);
	ip=*coor;jp=coor[1];kp=coor[2];

	// Set the global index ranges and size of this processor's domain
	ai=m*ip/mp;bi=m*(ip+1)/mp;sm=bi-ai;
	aj=n*jp/np;bj=n*(jp+1)/np;sn=bj-aj;
	ak=o*kp/op;bk=o*(kp+1)/op;so=bk-ak;
	smn=sm*sn;
	smno=smn*so;
}

/** Initializes an array of length 27 with the ranks of neighboring processors
 * in a 3x3x3 grid surrounding this processor. If a neighbor doesn't exist, the
 * array entry is set to -1.
 * \param[in] neigh a pointer to the array to initialize. */
void geometry::set_up_neighbors(int *neigh) {
	int o[3],&i=*o,&j=o[1],&k=o[2];
	for(k=kp-1;k<=kp+1;k++) for(j=jp-1;j<=jp+1;j++) for(i=ip-1;i<=ip+1;i++,neigh++) {
		if((z_prd||(k>=0&&k<op))&&
		   (y_prd||(j>=0&&j<np))&&
		   (x_prd||(i>=0&&i<mp))) MPI_Cart_rank(cart,o,neigh);
		else *neigh=-1;
	}
}
