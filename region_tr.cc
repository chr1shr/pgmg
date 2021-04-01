/** \file region_tr.hh
 * \brief Function implementations for the region class, representing part of a
 * particular grid level within a multigrid hierarchy. This file only contains routines
 * involved in transferring between one grid structure and another. */

#include "region.hh"
#include <cstring>
#include <cmath>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
	int dgetrs_(char *trans_,int *n,int *nrhs,double *a,int *lda,
		    int *ipiv,double *b,int *ldb,int *info);
	int dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
}

/** This constructor sets up a region as part of a multigrid hierarchy, for the
 * special case when this grid has the same mathematical structure as the
 * parent grid, but is redistributed onto a new (smaller) processor geometry.
 * \param[in] *p a pointer to the parent region class.
 * \param[in] gm a reference to the new geometry.
 * \param[in] com_ a reference to a communication buffer class. */
region::region(region* p_,geometry &gm,comm_buffer &com_)
	: p(p_), rank(gm.rank), x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
	gs_enable(true), m(gm.m), n(gm.n), o(gm.o),
	ip(gm.ip), jp(gm.jp), kp(gm.kp), mp(gm.mp), np(gm.np), op(gm.op),
	cart(gm.cart), pcart(p->cart), At(NULL), ai(gm.ai), aj(gm.aj), ak(gm.ak),
	bi(gm.bi), bj(gm.bj), bk(gm.bk), sm(gm.sm), sn(gm.sn), so(gm.so),
	smn(gm.smn), smno(gm.smno), c_inf(NULL), c_ptr(NULL), i_inf(NULL),
	i_ptr(NULL), tr_size(0), tr_inf(NULL), Aexact(NULL), com(com_) {

	grid_message(" [transfer]");

	// Set up neighbor table and initialize memory
	gm.set_up_neighbors(neighbor);
	setup_communication();
	S=new int[6*smno];
	A=new double*[smno];

	// Set up the incoming transfer table, and use to to receive the RAT
	// bounds
	setup_incoming_table();

	// Allocate memory for problem entries, set up pointers, and compute
	// ghost grid size
	int i,j,k,s0,s1,s2,s3,s4,s5,*Sp=S;
	double *Amp=(Am=new double[Asize]),**Ap=A;
	li=ai;hi=bi;lj=aj;hj=bj;lk=ak;hk=bk;
	for(k=ak;k<bk;k++) for(j=aj;j<bj;j++) for(i=ai;i<bi;i++) {
		s0=*(Sp++);if(i+s0<li) li=i+s0;
		s1=*(Sp++);if(i+s1>hi) hi=i+s1;
		s2=*(Sp++);if(j+s2<lj) lj=j+s2;
		s3=*(Sp++);if(j+s3>hj) hj=j+s3;
		s4=*(Sp++);if(k+s4<lk) lk=k+s4;
		s5=*(Sp++);if(k+s5>hk) hk=k+s5;
		*(Ap++)=Amp;
		Amp+=(s1-s0)*(s3-s2)*(s5-s4);
	}
	mg=hi-li;ng=hj-lj;og=hk-lk;mng=mg*ng;

	// Allocate function and source arrays and set up communication buffers
	x=new double[mng*og];
	r0=new double[smno];
	setup_gs_buffers();

	// Set up in the x and r pointers associated with the incoming transfer
	// table
	setup_incoming_pointers();

	// Check size for output matrix
	setup_output();
}

/** Sets up an incoming grid transfer table, for the case when this region is
 * part of a new grid geometry.
 * \param[in] p a pointer to the region in the parent grid, which will transfer
 *		information to this grid. */
void region::setup_incoming_table() {

	// Determine the range of processors in the new geometry that overlap
	// with this region's domain
	int ipl,ipu,jpl,jpu,kpl,kpu,
	    ui,uj,uk,vi,vj,vk,q[4],&i=*q,&j=q[1],&k=q[2],uijk,
	    *wsm=new int[p->mp+p->np+p->op],*wsn=wsm+p->mp,*wso=wsn+p->np;
	incoming_index_range(p->osm,wsm,p->mp,ai,bi,ipl,ipu);
	incoming_index_range(p->osn,wsn,p->np,aj,bj,jpl,jpu);
	incoming_index_range(p->oso,wso,p->op,ak,bk,kpl,kpu);

	// Allocate space for the transfer table and pointers
	allocate_transfer_memory((ipu-ipl)*(jpu-jpl)*(kpu-kpl));

	// Create the transfer table entries by looping over the range of
	// processors in the new geometry that overlap with this region's
	// domain
	double ***tAp=tr_A;
	int **tSp=tr_S,*tp=tr_inf;
	for(k=kpl;k<kpu;k++) {
		uk=k==kpl?0:wso[k-1]-ak;
		vk=k==kpu-1?so:wso[k]-ak;
		for(j=jpl;j<jpu;j++) {
			uj=j==jpl?0:wsn[j-1]-aj;
			vj=j==jpu-1?sn:wsn[j]-aj;
			for(i=ipl;i<ipu;i++,tp+=6) {
				ui=i==ipl?0:wsm[i-1]-ai;
				vi=i==ipu-1?sm:wsm[i]-ai;

				// Set the table entry
				MPI_Cart_rank(pcart,q,tp);
				tp[2]=vi-ui;
				tp[3]=vj-uj;
				tp[4]=vk-uk;
				tp[5]=tp[2]*tp[3]*tp[4];
#if PGMG_VERBOSE >=2
				printf("%d: ITAB <%d %d %d> (%d,%d,%d) [%d,%d,%d] {%d}\n",
				       rank,i,j,k,ui,uj,uk,tp[2],tp[3],tp[4],*tp);
#endif

				// Set the pointers
				*(tAp++)=A+(uijk=ui+sm*(uj+sn*uk));
				*(tSp++)=S+6*uijk;
			}
		}
	}

	// Delete the temporary memory used to store the coordinates of the
	// domains in the parent grid
	delete [] wsm;

	// Check for enough space for the x, r, and S transfers, both for the
	// parent's communication and for this region's communication
	tr_psmno=p->smno;
	com.check_buf_int(6*(tr_psmno+smno));
	com.check_buf(tr_psmno+smno);

	// Receive the RAT bounds and check that there will be enough buffer
	// space to send all of the matrix entries
	p->send_S();
	receive_S();
	p->send_wait();

	// Check for enough space for A transfer, both for the parent's
	// communication and this region's communication
	tr_pAsize=p->Asize;
	com.check_buf(tr_pAsize+Asize);
}

/** Calculates the range of processors in a particular dimension that overlap with
 * an given global index range.
 * \param[in] os an array of processor sizes in this dimension.
 * \param[in] ws an array in which to store the global indices of the processor
 *		 boundaries, computed by summing the processor sizes. The array
 *		 is only computed up to the u index, since that is all that
 *		 will be needed in the subsequent computations.
 * \param[in] vp the total number of processors in this dimension.
 * \param[in] (a,b) the global index range to consider.
 * \param[out] (l,u) the processor range the overlaps with the given global
 *		     index range. */
void region::incoming_index_range(int *os,int *ws,int vp,int a,int b,int &l,int &u) {
	*ws=*os;
	l=0;

	// Compare the lower index with the cumulative processor boundary
	// indices, until the former exceeds or equals the latter
	while(a>=ws[l]) {
		ws[l+1]=ws[l]+os[l+1];
		l++;
		if(l==vp) p_fatal_error("Error in transfer indexing (lower index)",1);
	}

	// Compare the upper index with the cumulative processor boundary
	// indices, until the former exceeds the latter
	u=l;
	while(b>ws[u]) {
		if(u+1==vp) p_fatal_error("Error in transfer indexing (upper index)",1);
		ws[u+1]=ws[u]+os[u+1];
		u++;
	}
	u++;
}

/** Sets up the pointers to the x and r arrays needed for the incoming transfer
 * table. These must be initialized later than the transfer table itself,
 * becuase the x and r arrays are not available when the transfer table is
 * initialized. */
void region::setup_incoming_pointers() {
	double **txp=tr_x,**trp=tr_r,***tAp=tr_A,***tAe=tAp+tr_size;
	int uijk,uj,uk;

	// Loop over the entries in the transfer table
	while(tAp<tAe) {

		// Convert the A transfer pointer into an index
		uijk=int(*(tAp++)-A);

		// Set up the x and r pointers, taking into account that the x
		// array is indexed differently and has ghost points
		uj=uijk/sm%sn;
		uk=uijk/smn;
		*(trp++)=r0+uijk;
		*(txp++)=x0+(uijk+(mg-sm)*uj+(mng-smn)*uk);
	}
}

/** Sets a table to transfer the contents of this region to other processors in
 * a new geometry.
 * \param[in] gm the new geometry to consider. */
void region::setup_outgoing_table(geometry &gm) {

	// Give an error if a transfer table is already set up. This would cause
	// a memory clash, and in any case, it would always be redundant for
	// this region to be involved in two transfers (incoming and outgoing).
	if(tr_inf!=NULL)
		p_fatal_error("Attempt to initialize second transfer table",1);

	// Set the flag to disable Gauss--Seidel sweeps on this level
	gs_enable=false;

	// Determine the range of processors in the new geometry that overlap
	// with this region's domain
	int ipl=(ai*gm.mp+gm.mp-1)/m,ipu=(bi*gm.mp-1)/m+1,
	    jpl=(aj*gm.np+gm.np-1)/n,jpu=(bj*gm.np-1)/n+1,
	    kpl=(ak*gm.op+gm.op-1)/o,kpu=(bk*gm.op-1)/o+1,
	    ui,uj,uk,vi,vj,vk,q[4],&i=*q,&j=q[1],&k=q[2],uijk;

	// Allocate space for the transfer table and pointers
	com.check_buf(Asize);
	allocate_transfer_memory((ipu-ipl)*(jpu-jpl)*(kpu-kpl));

	// Create the transfer table entries by looping over the range of
	// processors in the new geometry that overlap with this region's
	// domain
	double **txp=tr_x,**trp=tr_r,***tAp=tr_A;
	int **tSp=tr_S,*tp=tr_inf;
	for(k=kpl;k<kpu;k++) {
		uk=k==kpl?0:k*gm.o/gm.op-ak;
		vk=k==kpu-1?so:(k+1)*gm.o/gm.op-ak;
		for(j=jpl;j<jpu;j++) {
			uj=j==jpl?0:j*gm.n/gm.np-aj;
			vj=j==jpu-1?sn:(j+1)*gm.n/gm.np-aj;
			for(i=ipl;i<ipu;i++,tp+=6) {
				ui=i==ipl?0:i*gm.m/gm.mp-ai;
				vi=i==ipu-1?sm:(i+1)*gm.m/gm.mp-ai;

				// Set the table entry
				*tp=gm.ranks[i+gm.mp*(j+gm.np*k)];
				tp[2]=vi-ui;
				tp[3]=vj-uj;
				tp[4]=vk-uk;
				tp[5]=tp[2]*tp[3]*tp[4];
#if PGMG_VERBOSE >= 2
				printf("%d: OTAB <%d %d %d> (%d,%d,%d) [%d,%d,%d] {%d}\n",
				       rank,i,j,k,ui,uj,uk,tp[2],tp[3],tp[4],*tp);
#endif

				// Set the pointers
				*(txp++)=x0+(ui+mg*(uj+ng*uk));
				*(trp++)=r0+(uijk=ui+sm*(uj+sn*uk));
				*(tAp++)=A+uijk;
				*(tSp++)=S+6*uijk;
			}
		}
	}
}

/** Allocates the memory for performing grid transfers.
 * \param[in] tr_size_ the number of transfers to different processors. */
void region::allocate_transfer_memory(int tr_size_) {

	// Allocate buffers for the transfer information and for the pointers
	// to the corner gridpoint in each transfer
	tr_size=tr_size_;
	tr_inf=new int[6*tr_size];
	tr_x=new double*[tr_size];
	tr_r=new double*[tr_size];
	tr_A=new double**[tr_size];
	tr_S=new int*[tr_size];

	// If the number of transfers exceeds the request/status arrays, then
	// allocate more space
	if(tr_size>req_stat_size) {
		delete [] req;
		delete [] stat;
		req_stat_size=tr_size;
		req=new MPI_Request[req_stat_size];
		stat=new MPI_Status[req_stat_size];
	}

	// Check buffer size is big enough for transfering the fields and the
	// RAT bounds
	com.check_buf_int(6*smno);
	com.check_buf(smno);
}

/** Initiates non-blocking sends of the source term of this region to the
 * relevant processors in the child region. */
void region::send_r() {
	double *pp=com.buf,*outp,**trp=tr_r;
	int i,j,k;
	MPI_Request *reqp=req;

	// Assemble the buffers and send them
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,reqp++,trp++) {
		outp=pp;
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++)
			*(pp++)=(*trp)[i+sm*(j+sn*k)];
		MPI_Isend(outp,tp[5],MPI_DOUBLE,*tp,msg_trans_r,cart,reqp);
	}
}

/** Initiates non-blocking sends of the solution on this region to the relevant
 * processors in the parent region. */
void region::send_x() {
	double *pp=com.buf+tr_psmno,*outp,**txp=tr_x;
	int i,j,k;
	MPI_Request *reqp=req;

	// Assemble the buffers and send them. Note that the parent region's
	// communicator is used, since this region's communicator may not
	// include all relevant processes.
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,reqp++,txp++) {
		outp=pp;
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++)
			*(pp++)=(*txp)[i+mg*(j+ng*k)];
		MPI_Isend(outp,tp[5],MPI_DOUBLE,*tp,msg_trans_x,pcart,reqp);
	}
}

/** Receives the source term contributions from the parent regions. */
void region::receive_r() {
	int i,j,k;
	double *pp=com.buf+tr_psmno,**trp=tr_r;
	MPI_Request *reqp=req;

	// Receive the source term contributions into the communication buffer.
	// Note that the parent region's communicator is used, since this
	// region's communicator may not include all relevant processes.
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;reqp++,pp+=tp[5],tp+=6)
		MPI_Irecv(pp,tp[5],MPI_DOUBLE,*tp,msg_trans_r,pcart,reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the source term array
	pp=com.buf+tr_psmno;
	MPI_Waitall(tr_size,req,stat);
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,trp++)
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++)
			(*trp)[i+sm*(j+sn*k)]=*(pp++);
}

/** Receives the solution from the child regions. */
void region::receive_x() {
	int i,j,k;
	double *pp=com.buf,**txp=tr_x;
	MPI_Request *reqp=req;

	// Receive the solution contributions from the child regions into the
	// communication buffer
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;reqp++,pp+=tp[5],tp+=6)
		MPI_Irecv(pp,tp[5],MPI_DOUBLE,*tp,msg_trans_x,cart,reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the solution array
	pp=com.buf;
	MPI_Waitall(tr_size,req,stat);
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,txp++)
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++)
			(*txp)[i+mg*(j+ng*k)]=*(pp++);
}

/** Initiates non-blocking sends of the RAT bound information to the child
 * regions as part of a grid transfer. The routine also scans the RAT bound
 * information to calculate the how many RAT entries will be subsequently
 * communicated. */
void region::send_S() {
	int *pp=reinterpret_cast<int*>(com.buf),*outp,**tSp=tr_S,*Sp,*Se;
	int i,j,k;
	MPI_Request *reqp=req;

	// Assemble the buffers and send them, calculating the size of the
	// subsequent RAT entry messages and storing them in tp[1] in the
	// transfer information
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,reqp++,tSp++) {
		outp=pp;
		tp[1]=0;
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++) {
			Sp=*tSp+6*(i+sm*(j+sn*k));Se=Sp+6;
			tp[1]+=rat_size(Sp);
			while(Sp<Se) *(pp++)=*(Sp++);
		}
		MPI_Isend(outp,6*tp[5],MPI_INT,*tp,msg_trans_S,cart,reqp);
	}
}

/** Receives the RAT bound information from the parent regions. The routine
 * also scans the RAT bound information to calculate how many RAT entries will
 * be subsequently communicated. */
void region::receive_S() {
	int *ps=reinterpret_cast<int*>(com.buf)+6*tr_psmno,*pp=ps,**tSp=tr_S,*Sp,*Se;
	int i,j,k;
	MPI_Request *reqp=req;

	// Receive the RAT bound information into the communication buffer.
	// Note that the parent region's communicator is used, since this
	// region's communicator may not include all relevant processes.
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;reqp++,pp+=6*tp[5],tp+=6)
		MPI_Irecv(pp,6*tp[5],MPI_INT,*tp,msg_trans_S,pcart,reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the RAT bound array. In addition, calculate the size of the
	// subsequent RAT entry messages and store them in tp[1] in the
	// transfer information
	MPI_Waitall(tr_size,req,stat);
	pp=ps;Asize=0;
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,tSp++) {
		tp[1]=0;
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++) {
			Sp=*tSp+6*(i+sm*(j+sn*k));Se=Sp+6;
			tp[1]+=rat_size(pp);
			while(Sp<Se) *(Sp++)=*(pp++);
		}
		Asize+=tp[1];
	}
}

/** Initiates non-blocking sends of the RAT matrix entries to the child regions,
 * as part of a grid transfer. */
void region::send_A() {
	double *pp=com.buf,*outp,***tAp=tr_A,*Ap,*Ae;
	int **tSp=tr_S,i,j,k,ijk;
	MPI_Request *reqp=req;

	// Assemble the buffers and send them
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,reqp++,tSp++,tAp++) {
		outp=pp;
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++) {
			Ap=(*tAp)[ijk=i+sm*(j+sn*k)];
			Ae=Ap+rat_size(*tSp+6*ijk);
			while(Ap<Ae) *(pp++)=*(Ap++);
		}
		MPI_Isend(outp,tp[1],MPI_DOUBLE,*tp,msg_trans_A,cart,reqp);
	}
}

/** Receives the RAT matrix entries from the parent regions, as part of a grid
 * transfer. */
void region::receive_A() {
	double *pp=com.buf+tr_pAsize,***tAp=tr_A,*Ap,*Ae;
	int **tSp=tr_S,i,j,k,ijk;
	MPI_Request *reqp=req;

	// Receive the RAT matrix entries into the communication buffer. Note
	// that the parent region's communicator is used, since this region's
	// communicator may not include all relevant processes.
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;reqp++,pp+=tp[1],tp+=6)
		MPI_Irecv(pp,tp[1],MPI_DOUBLE,*tp,msg_trans_A,pcart,reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the RAT matrix entry array
	pp=com.buf+tr_pAsize;
	MPI_Waitall(tr_size,req,stat);
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,tSp++,tAp++) {
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++) {
			Ap=(*tAp)[ijk=i+sm*(j+sn*k)];
			Ae=Ap+rat_size(*tSp+6*ijk);
			while(Ap<Ae) *(Ap++)=*(pp++);
		}
	}
}

/** Attempts to switch on exact computation for this region. This will only
 * occur if this process has no neighbors, and if the total problem size is
 * small enough. */
void region::enable_exact() {
	if(mp==1&&np==1&&op==1&&smno<pgmg_max_exact) {
		ipiv=new int[smno];
		Aexact=new double[smno*smno];
		com.check_buf(smno);
	}
}

/** Assuming exact computation is enabled for this region, this routine
 * performs the LU decomposition of the linear system, allowing for it to be
 * subsequently solved exactly. */
void region::lu_decomposition() {

	// Clear the exact linear system array
	double *pp=Aexact,*pe=pp+smno*smno,*Amp,**Ap=A;
	while(pp<pe) *(pp++)=0;

	// Construct the dense matrix to solve
	pp=Aexact;
	int i,j,k,di,dj,dk,ei,ej,ek,*Sp=S;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++,Sp+=6,pp++) {
		Amp=*(Ap++);
		for(dk=k+Sp[4];dk<k+Sp[5];dk++) {
			ek=step_mod(dk,so);
			for(dj=j+Sp[2];dj<j+Sp[3];dj++) {
				ej=step_mod(dj,sn);
				for(di=i+*Sp;di<i+Sp[1];di++) {
					ei=step_mod(di,sm);
					pp[((ek*sn+ej)*sm+ei)*smno]=*(Amp++);
				}
			}
		}
	}
    
    /*
    if(rank==0) for(i=0;i<smno;i++) {
        printf("%d [", i);
        for(j=0;j<smno;j++) {
            printf("%g ", Aexact[j + i*smno]);
        }
        printf("]\n");
    }
    */

	// Perform the LU decomposition using LAPACK
	int info;
	dgetrf_(&smno,&smno,Aexact,&smno,ipiv,&info);
//	if (info!=0) p_fatal_error("LAPACK LU decomposition failed",1);

//	if (rank == 0 && info == 0) 
//		for (int i = 0; i < smno; i++) printf("%.12g\n",Aexact[i*(1+smno)]);

	// Check the conditioning of the linear system. If it's poorly
	// conditioned, then just fall back on Gauss--Seidel.
	double eig_hi=fabs(*Aexact),eig_lo=fabs(Aexact[smno*smno-1]);
	if (info != 0 || eig_hi==0 || eig_lo/eig_hi<pgmg_exact_tol) {
		if (rank == 0) printf("#-->LAPACK LU decomposition failed\n"
			"#   Falling back on Gauss-Seidel\n");
		delete [] Aexact;
		delete [] ipiv;
		Aexact=NULL;
	}
}

/** Attempts to solve the linear system exactly, assuming that an LU
 * decomposition of the system has already been set up. If no LU decomposition
 * is available, the routine falls back on performing a large number of
 * Gauss--Seidel sweeps. */
void region::solve_exact(bool smooth) {
	if(Aexact!=NULL && !smooth) {

		// Copy the source term into temporary work space
		memcpy(com.buf,r0,smno*sizeof(double));

		// Solve the linear system using the previously computed LU
		// decomposition
		char trans='N';
		int info,nrhs=1,i,j,k;
		dgetrs_(&trans,&smno,&nrhs,Aexact,&smno,ipiv,com.buf,&smno,&info);

		// Copy the result of the linear solve into the main solution array
		double *pp=com.buf;
		for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++){
 			x0[(k*ng+j)*mg+i]=*(pp++);}

	} else for(int l=0;l<pgmg_gs_exact;l++) gauss_seidel();
}

/** Prints some statistics about the sums of the fields. */
void region::diagnostic_sums() {
	double su[12];

	// Make sure the ghost regions are up to date
	communicate();

	// Clear the accumulators
	int i,j,k;
	for(i=0;i<6;i++) su[i]=0.;

	// Loop over the primary points in the grid
	double q;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		*su+=(q=residual(i,j,k));
		su[1]+=q*q;
		su[2]+=(q=r0[(k*sn+j)*sm+i]);
		su[3]+=q*q;
		su[4]+=(q=x0[(k*ng+j)*mg+i]);
		su[5]+=q*q;
	}

	MPI_Reduce(su,su+6,6,MPI_DOUBLE,MPI_SUM,0,cart);
	if(rank==0) printf("Grid (%d,%d,%d): res (%g %g) r(%g %g) x(%g %g)\n",
			   m,n,o,su[6],su[7],su[8],su[9],su[10],su[11]);
}
