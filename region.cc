/** \file region.hh
 * \brief Function implementations for the region class, representing part of a
 * particular grid level within a multigrid hierarchy. */

#include "region.hh"
#include <cstring>

/** This constructor sets up a region for a lower level of a multigrid computation,
 * by linking it to the parent region class.
 * \param[in] *p_ a pointer to the parent region class.
 * \param[in] com_ a reference to a communication buffer class. */
region::region(region *p_,comm_buffer &com_)
	: p(p_), rank(p->rank), x_prd(p->x_prd), y_prd(p->y_prd), z_prd(p->z_prd),
	gs_enable(true), m((p->m+(x_prd?1:2))>>1), n((p->n+(y_prd?1:2))>>1), o((p->o+(z_prd?1:2))>>1),
	ip(p->ip), jp(p->jp), kp(p->kp), mp(p->mp), np(p->np), op(p->op),
	cart(p_->cart), pcart(MPI_COMM_NULL), At(new double[init_A_temp_size]),
	Atmem(init_A_temp_size), ai((p->ai+1)>>1), aj((p->aj+1)>>1),
	ak((p->ak+1)>>1), bi(ip==mp-1?m:(p->bi+1)>>1),
	bj(jp==np-1?n:(p->bj+1)>>1), bk(kp==op-1?o:(p->bk+1)>>1),
	sm(bi-ai), sn(bj-aj), so(bk-ak), smn(sm*sn), smno(smn*so),
	pmg(p->mg), png(p->ng), pmng(pmg*png),
	px0(p->x0+(p->ai&1)+(p->aj&1?pmg:0)+(p->ak&1?pmng:0)),
	pS0(p->S+6*((p->ai&1)+(p->aj&1?(p->sm):0)+(p->ak&1?p->sm*p->sn:0))),
	c_inf(NULL), c_ptr(NULL), i_inf(NULL), i_ptr(NULL), tr_size(0), tr_inf(NULL),
	Aexact(NULL), com(com_), out_of_bounds(-99) {

	grid_message("");

	// Set up neighbor table and initialize memory
	memcpy(neighbor,p->neighbor,27*sizeof(int));
	setup_communication();
	S=new int[6*smno];
	A=new double*[smno];

	// Set up the strip buffers
	setup_strip_buffers();

	// Compute RAT matrix bounds using the parent region
	fill_rat_bounds();

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

	// Set up the pointers used to communicate the interpolation and
	// restriction strips. Since they point at the solution and source term
	// arrays, they have to be initialized here, after the arrays has been
	// allocated.
	setup_rt_pointers();

	// Check size for output matrix
	setup_output();
}

/** The class destructor deallocates the dynamically allocated arrays. */
region::~region() {

	// Delete the exact solution arrays if present
	if(Aexact!=NULL) {
		delete [] Aexact;
		delete [] ipiv;
	}

	// Delete grid transfer arrays if present
	if(tr_inf!=NULL) {
		delete [] tr_S;delete [] tr_A;
		delete [] tr_r;delete [] tr_x;
		delete [] tr_inf;
	}

	// Delete output information array
	delete [] osm;

	// Delete communication arrays if present
	if(c_inf!=NULL) {
		delete [] c_inf;delete [] c_ptr;
	}

	// Delete restriction and interpolation arrays if present
	if(i_inf!=NULL) {
		delete [] Atrans;delete [] Strans;
		delete [] S2_ptr;delete [] i2_ptr;
		delete [] i2_inf;delete [] Ac_size;
		delete [] A_ptr;delete [] S_ptr;
		delete [] r_ptr;delete [] i_ptr;
		delete [] i_inf;
	}

	// Delete primary arrays and MPI status/request arrays
	delete [] r0;delete [] x;
	delete [] Am;delete [] S;delete [] A;
	delete [] stat;delete [] req;

	// Delete the local RAT computation array
	if(At!=NULL) delete [] At;
}

int region::gs_num = 0;

/** Sets up the range of neighbors to communicate with, and the MPI request and
 * status arrays. */
void region::setup_communication() {

	// Set up block of neighbors to contact
	lkp=neighbor[4]==-1?1:0;
	ljp=neighbor[10]==-1?1:0;
	lip=neighbor[12]==-1?1:0;
	hip=neighbor[14]==-1?2:3;
	hjp=neighbor[16]==-1?2:3;
	hkp=neighbor[22]==-1?2:3;

	// Compute the total number of neighbors
	tneigh=(hip-lip)*(hjp-ljp)*(hkp-lkp)-1;
	// Set up the MPI request and status arrays
	req_stat_size=max(rank==0?mp*np*op+1:1,2*tneigh);
	req=new MPI_Request[req_stat_size];
	stat=new MPI_Status[req_stat_size];
}

/** Sets up the buffers used in the Gauss-Seidel sweep, and communicates
 * the buffer sizes between neighboring processes. */
void region::setup_gs_buffers() {

	// If there are no neighbors, then just return
	if(tneigh==0) {
		cneigh_in=cneigh_out=0;
		x0=x;
		indent = 0;
		return;
	}

	// Set up pointer to (0,0,0) element of the x array
	int xd_size=ai-li,yd_size=aj-lj,zd_size=ak-lk;
	indent = (xd_size+mg*(yd_size+ng*zd_size));
	x0=x+indent;

	// Allocate memory and transfer buffer sizes to neighbors
	int *ctmp=new int[12*tneigh];
	transfer_buffer_sizes(ctmp,xd_size,hi-bi,yd_size,hj-bj,zd_size,hk-bk,cneigh_in,cneigh_out);

	// Set up the memory for the non-zero buffers
	c_inf=new int[6*(cneigh_in+cneigh_out)];
	c_ptr=new double*[cneigh_in+cneigh_out];

	// Copy information about incoming buffers to new data structure
	int *cp=ctmp,*cp2=c_inf,i,j,k;
	double **pp=c_ptr;
	for(k=hkp-1;k>=lkp;k--) for(j=hjp-1;j>=ljp;j--) for(i=hip-1;i>=lip;i--) {
		if(i==1&&j==1&&k==1) continue;
		if(cp[5]>0) {
			*(pp++)=x0+(i==1?0:(i!=2?-xd_size:sm))
				  +mg*((j==1?0:(j!=2?-yd_size:sn))
				  +ng*(k==1?0:(k!=2?-zd_size:so)));
			six_copy(cp2,cp);
		} else cp+=6;
	}

	// Copy information about outgoing buffers to new data structure
	for(k=lkp;k<hkp;k++) for(j=ljp;j<hjp;j++) for(i=lip;i<hip;i++) {
		if(i==1&&j==1&&k==1) continue;
		if(cp[5]>0) {
			*(pp++)=x0+(i!=2?0:sm-cp[2])
				  +mg*((j!=2?0:sn-cp[3])
				  +ng*(k!=2?0:so-cp[4]));
			six_copy(cp2,cp);
		} else cp+=6;
	}

	// Delete temporary buffer
	delete [] ctmp;
}

/** Exchanges the sizes of the ghost regions for a particular type of
 * communication with the neighboring processors, and creates information about
 * the messages that need to be sent and received, eliminating any messages
 * that would have zero size.
 * \param[in] ctmp an array in which to assemble the information about the
 *		   messages.
 * \param[in] (xd_size,xu_size) the region widths in the negative and positive
 *				x directions, respectively.
 * \param[in] (yd_size,yu_size) the region widths in the negative and positive
 *				y directions, respectively.
 * \param[in] (zd_size,zu_size) the region widths in the negative and positive
 *				z directions, respectively.
 * \param[out] n_in the total number of messages that will be received from
 *		    other processors.
 * \param[out] n_out the total number of message that will be sent to other
 *		     processors. */
int region::transfer_buffer_sizes(int *ctmp,int xd_size,int xu_size,int yd_size,int yu_size,int zd_size,int zu_size,int &n_in,int &n_out) {
	int i,j,k,ijk,l=0,*cp=ctmp;
	MPI_Request *reqp=req;

	for(k=hkp-1;k>=lkp;k--) for(j=hjp-1;j>=ljp;j--) for(i=hip-1;i>=lip;i--) {
		if((ijk=i+3*(j+3*k))==13) continue;

		// Store neighbor information at the size of the communication box
		*(cp++)=neighbor[ijk];
		*(cp++)=ijk;
		*cp=i==1?sm:(i!=2?xd_size:xu_size);
		cp[1]=j==1?sn:(j!=2?yd_size:yu_size);
		cp[2]=k==1?so:(k!=2?zd_size:zu_size);
		l+=cp[3]=*cp*cp[1]*cp[2];

		// Store a pointer to the box corner
		MPI_Isend(cp,4,MPI_INT,neighbor[ijk],msg_tbs|ijk,cart,reqp);
		cp+=4;reqp++;
	}

	// Receive the output buffer sizes from the neighbors
	for(k=lkp;k<hkp;k++) for(j=ljp;j<hjp;j++) for(i=lip;i<hip;i++) {
		if((ijk=i+3*(j+3*k))==13) continue;

		*(cp++)=neighbor[ijk];
		*(cp++)=26-ijk;
		MPI_Irecv(cp,4,MPI_INT,neighbor[ijk],msg_tbs|(26-ijk),cart,reqp);

		// Store a pointer to the box corner
		cp+=4;reqp++;
	}

	// Wait for all sends to complete
	MPI_Waitall(2*tneigh,req,stat);

	// Count non-zero communication buffers, and complete the total
	// communication memory count
	n_in=n_out=0;
	cp=ctmp+5;
	for(;cp<ctmp+6*tneigh;cp+=6) if(*cp>0) n_in++;
	for(;cp<ctmp+12*tneigh;cp+=6) {
		l+=*cp;
		if(*cp>0) n_out++;
	}

	// Check that there is enough memory in the communication buffer for
	// all incoming and outgoing region communications
	com.check_buf(l);
	return l;
}

/** Sets up the strip buffers that are used in interpolation, restriction,
 * and RAT computation. */
void region::setup_strip_buffers() {

	// Calculate the widths of the ghost strips
	int xd_size=(p->ai)&1,xu_size=neighbor[14]==-1?0:((p->bi&1)==0?1:0),
	    yd_size=(p->aj)&1,yu_size=neighbor[16]==-1?0:((p->bj&1)==0?1:0),
	    zd_size=(p->ak)&1,zu_size=neighbor[22]==-1?0:((p->bk&1)==0?1:0),q;

	// Communicate the widths to the neighbors, storing the results in
	// a temporary buffer
	int *ctmp=new int[12*tneigh];
	ntrans=transfer_buffer_sizes(ctmp,xd_size,xu_size,yd_size,yu_size,
				     zd_size,zu_size,ineigh_in,ineigh_out);

	// Set up the memory for the non-zero buffers
	i_inf=new int[6*(ineigh_in+ineigh_out)];
	i_ptr=new double*[ineigh_out];
	r_ptr=new double*[ineigh_out];
	S_ptr=new int*[ineigh_out];
	A_ptr=new double**[ineigh_out];
	Ac_size=new int[ineigh_in+ineigh_out];

	// Set up memory for the buffer look-up tables
	i2_inf=new int[27*4];
	i2_ptr=new double*[27];
	S2_ptr=new int*[27];
	for(q=3;q<27*4;q+=4) i2_inf[q]=-1;

	// Set up memory to save the RAT contribution dimensions that are
	// communicated to neighbors. These need to be stored in order to send
	// the RAT contributions at a later stage.
	Strans=new int[6*ntrans];
	Atrans=new double*[ntrans];

	// Copy information about incoming buffers to new data structure
	int *cp=ctmp,*cp2=i_inf,i,j,k;
	for(k=hkp-1;k>=lkp;k--) for(j=hjp-1;j>=ljp;j--) for(i=hip-1;i>=lip;i--) {
		if(i==1&&j==1&&k==1) continue;
		if(cp[5]>0) {

			// Set up the interpolation look-up table
			int *i2p=i2_inf+4*(i+3*j+9*k);
			*(i2p++)=cp[2];
			*(i2p++)=cp[3];
			*(i2p++)=(i==1?0:(i!=2?-xd_size:sm))
				 +cp[2]*((j==1?0:(j!=2?-yd_size:sn))
				 +cp[3]*(k==1?0:(k!=2?-zd_size:so)));
			*i2p=cp[5];
			six_copy(cp2,cp);
		} else cp+=6;
	}

	// Copy information about outgoing buffers to new data structure
	int **pp=S_ptr,ijk;
	double ***p2=A_ptr;
	for(k=lkp;k<hkp;k++) for(j=ljp;j<hjp;j++) for(i=lip;i<hip;i++) {
		if(i==1&&j==1&&k==1) continue;
		if(cp[5]>0) {

				// Set pointers for receiving RAT strip contributions
			ijk=(i!=2?0:sm-cp[2])+sm*((j!=2?0:sn-cp[3])+sn*(k!=2?0:so-cp[4]));
			*(pp++)=S+6*ijk;
			*(p2++)=A+ijk;

			// Copy data into new data structure
			six_copy(cp2,cp);
		} else cp+=6;
	}

	// Delete temporary buffer
	delete [] ctmp;
}

/** Sets up pointers to the solution array that are filled by the interpolation
 * strip communication. These must be set up here, since the solution array is
 * not allocated when the strips are first determined. */
void region::setup_rt_pointers() {
	double **pp=i_ptr,**p2=r_ptr;
	for(int *cp=i_inf+6*ineigh_in,*ce=cp+6*ineigh_out;cp<ce;cp+=6) {
		// Use message tag to determine which direction this messages
		// is being sent to
		int i=cp[1]%3,j=(cp[1]/3)%3,k=cp[1]/9;

		// Set the pointers to the lowest corner of the grid points to
		// send
		*(pp++)=x0+(i!=0?0:sm-cp[2])
			+mg*((j!=0?0:sn-cp[3])
			+ng*(k!=0?0:so-cp[4]));
		*(p2++)=r0+(i!=0?0:sm-cp[2])
			+sm*((j!=0?0:sn-cp[3])
			+sn*(k!=0?0:so-cp[4]));
	}
}

/** Performs a Gauss--Seidel sweep, communicating edges to neighboring
 * processors as needed. */
void region::gauss_seidel(int mode) {
    gs_num ++;
	int i,j,k,*Sp=S,is,js,id,jd;
	double Ax,Axc,*Amp,**Ap=A,*xcp,*xp,*xkp,*xjp,*xip,*rp=r0;

	// Fill in ghost cells with values from neighboring processors
	communicate();

#if defined(DEBUG)
    bool hasnan=false;
#endif
	// Carry out forward Gauss--Seidel sweep
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {

		Ax=0;
		// Compute pointers and jumps
		Amp=*(Ap++);
		xcp=x0+(i+mg*(j+ng*k));
		xkp=xcp+*Sp+mg*Sp[2];xp=xkp+Sp[4]*mng;
		is=Sp[1]-*Sp;js=(Sp[3]-Sp[2])*mg;
		id=mg-is;jd=mng-js;
		// Compute the terms up to the central element
		while(xp!=xkp) {
			xjp=xp+js;
			while(xp!=xjp) {
				xip=xp+is;
				while(xp!=xip) {
					Ax+=(*(Amp++))*(*(xp++));
				}
				xp+=id;
			}
			xp+=jd;
		}
		xkp=xp+Sp[5]*mng;
		xjp=xp-Sp[2]*mg;
		while(xp!=xjp) {
			xip=xp+is;
			while(xp!=xip){
				Ax+=(*(Amp++))*(*(xp++));
			}
			xp+=id;
		}
		xjp=xp+Sp[3]*mg;
		while(xp!=xcp){
			Ax+=(*(Amp++))*(*(xp++));
		}
		xip=xp+Sp[1];

		// Remember the central element
		Axc=1./(*(Amp++));xp++;
#if defined(DEBUG)
        if(std::isnan(Axc)) {
            printf("Region:: shit central element is zero\n");
            hasnan=true;
        }
#endif

		// Compute the rest of the terms
		while(xp!=xip) {
			Ax+=(*(Amp++))*(*(xp++));
		}
		xp+=id;
		while(xp!=xjp) {
			xip=xp+is;
			while(xp!=xip){
				Ax+=(*(Amp++))*(*(xp++));
			}
			xp+=id;
		}
		xp+=jd;
		while(xp!=xkp) {
			xjp=xp+js;
			while(xp!=xjp) {
				xip=xp+is;
				while(xp!=xip){
					Ax+=(*(Amp++))*(*(xp++));
				}
				xp+=id;
			}
			xp+=jd;
		}

		// Perform the Gauss--Seidel update
		*xcp=(*(rp++)-Ax)*Axc;
		Sp+=6;
#if defined(DEBUG)
        if(std::isnan(*xcp)) {
            //printf("Region:: guass seidel result at (%d %d %d) is nan\n", i, j, k);
            hasnan=true;
        }
#endif
	}

#if defined(DEBUG)
    if(hasnan) p_fatal_error("Region:: nans encountered in Gauss-Seidel",1);
#endif

	// If we are doing regular Gauss-Seidel, return
	if(!mode) return;
	// Otherwise, assume we are doing symmetric Gauss-Seidel

	// Fill in ghost cells with values from neighboring processors
	communicate();

	// Note that the ptrs to S, A are all the way at the back
	// We just need to rewind before using them

	// Carry out forward Gauss-Seidel sweep
	for(k=so-1;k>=0;k--) for(j=sn-1;j>=0;j--) for(i=sm-1;i>=0;i--) {

		Ax=0;
		// Compute pointers and jumps
		Amp=*(--Ap);
		Sp-=6;
		// Central element
		xcp=x0+(i+mg*(j+ng*k));
		xkp=xcp+*Sp+mg*Sp[2];xp=xkp+Sp[4]*mng;
		is=Sp[1]-*Sp;js=(Sp[3]-Sp[2])*mg;
		id=mg-is;jd=mng-js;
		// Compute the terms up to the central element
		while(xp!=xkp) {
			xjp=xp+js;
			while(xp!=xjp) {
				xip=xp+is;
				while(xp!=xip) {
					Ax+=(*(Amp++))*(*(xp++));
				}
				xp+=id;
			}
			xp+=jd;
		}
		xkp=xp+Sp[5]*mng;
		xjp=xp-Sp[2]*mg;
		while(xp!=xjp) {
			xip=xp+is;
			while(xp!=xip){
				Ax+=(*(Amp++))*(*(xp++));
			}
			xp+=id;
		}
		xjp=xp+Sp[3]*mg;
		while(xp!=xcp){
			Ax+=(*(Amp++))*(*(xp++));
		}
		xip=xp+Sp[1];

		// Remember the central element
		Axc=1./(*(Amp++));xp++;

		// Compute the rest of the terms
		while(xp!=xip) {
			Ax+=(*(Amp++))*(*(xp++));
		}
		xp+=id;
		while(xp!=xjp) {
			xip=xp+is;
			while(xp!=xip){
				Ax+=(*(Amp++))*(*(xp++));
			}
			xp+=id;
		}
		xp+=jd;
		while(xp!=xkp) {
			xjp=xp+js;
			while(xp!=xjp) {
				xip=xp+is;
				while(xp!=xip){
					Ax+=(*(Amp++))*(*(xp++));
				}
				xp+=id;
			}
			xp+=jd;
		}

		// Perform the Gauss--Seidel update
		*xcp=(*(--rp)-Ax)*Axc;
#if defined(DEBUG)
        if(std::isnan(*xcp)) {
            //printf("Region:: guass seidel result at (%d %d %d) is nan\n", i, j, k);
            hasnan=true;
        }
#endif
	}

#if defined(DEBUG)
    if(hasnan) p_fatal_error("Region:: nans encountered in Gauss-Seidel",1);
#endif
}

/** Performs a restriction operation, by coarsening the residual on the parent
 * grid in the multigrid hierarchy. */
void region::restriction() {

	// Obtain the ghost points from the neighboring processors, since they
	// will be needed to construct the residuals to be restricted
	p->communicate();

	// Set pointers to outgoing restriction buffers and clear them
	double *pp=com.buf;
	for(int q=26;q>=0;q--)
		if(i2_inf[4*q+3]>0) {i2_ptr[q]=pp;pp+=i2_inf[4*q+3];} else i2_ptr[q]=NULL;
	for(double *p2=com.buf;p2<pp;p2++) *p2=0.;

	// Set the central element of the interpolation information to point at
	// the source term
	i2_inf[52]=sm;i2_inf[53]=sn;i2_inf[54]=0;i2_ptr[13]=r0;

	// Clear the source term grid. TODO: it would be nice not to have to do
	// this as a separate step, but it may cause problems for the rest of
	// the routine.
	for(double *rp=r0,*re=r0+smno;rp<re;rp++) *rp=0.;

	// Scan the grid and fill in the restriction values
	rt_scan_grid(1024);

	// Send the restriction
	communicate_restriction_strips();
}

/** Performs an interpolation operation. The solution on this grid is
 * trilinearly interpolated and added to the solution on the parent grid. */
void region::interpolation() {

	// Communicate the edge strips of the vector field that are needed in
	// the interpolation
	communicate_interpolation_strips();

	// Set pointers to interpolation buffers
	double *pp=com.buf;
	for(int q=26;q>=0;q--)
		if(i2_inf[4*q+3]>0) {i2_ptr[q]=pp;pp+=i2_inf[4*q+3];} else i2_ptr[q]=NULL;

	// Set the central element of the interpolation information to point at
	// the solution vector
	i2_inf[52]=mg;i2_inf[53]=ng;i2_inf[54]=0;i2_ptr[13]=x0;

	// Scan the grid and fill in the interpolation values
	rt_scan_grid(0);
}

/** Calculates the bounds of the linear system at this level using the RAT
 * calculation, communicating with neighboring processors to obtain required
 * information about border cases. */
void region::fill_rat_bounds() {

	// Set pointers to outgoing RAT buffers and clear them
	int *pp=Strans;
	for(int q=26;q>=0;q--)
		if(i2_inf[4*q+3]>0) {S2_ptr[q]=pp;pp+=6*i2_inf[4*q+3];} else S2_ptr[q]=NULL;
	for(int *p2=Strans;p2<pp;) {*(p2++)=bound_high;*(p2++)=bound_low;}

	// Set the central element of the communication pointer to the main
	// bound array
	i2_inf[52]=sm;i2_inf[53]=sn;i2_inf[54]=0;S2_ptr[13]=S;

	// Clear the bounds
	for(int *Sp=S,*Se=S+6*smno;Sp<Se;) {*(Sp++)=0;*(Sp++)=1;}

	// Scan the grid and fill in the restriction values
	rt_scan_grid(3072);

	// Send the contributions to other processors
	communicate_rat_bound_strips();

	// Count the total amount of memory needed for the linear system
	// coefficients and return it
	Asize=rat_size(S);
	for(int *Sp=S+6,*Se=S+6*smno;Sp<Se;Sp+=6) Asize+=rat_size(Sp);
}

/** Calculates the linear system coefficients at the level using the RAT
 * calculation, communicating with neighboring processors to obtain required
 * information about border cases. It assumes that the sizes of the linear
 * system have already been established using the fill_rat_bounds routine. */
void region::compute_rats() {

	// Set pointers to outgoing RAT buffers and clear them. Note that
	// this is a tiny bit wasteful, since the coefficient storage is
	// recomputed from the Strans information.
	double *pp=com.buf,**p2=Atrans;
	int *Sp=Strans,*Se;
	for(int q=26;q>=0;q--) {
		if(i2_inf[4*q+3]>0) {
			i2_ptr[q]=reinterpret_cast<double*>(p2);
			for(Se=Sp+6*i2_inf[4*q+3];Sp<Se;*(p2++)=pp,pp+=rat_size(Sp),Sp+=6);
		} else i2_ptr[q]=NULL;
	}
	for(double *p3=com.buf;p3<pp;p3++) *p3=0.;

	// Set the central element of the communication pointer to the main
	// bound array
	i2_inf[52]=sm;i2_inf[53]=sn;i2_inf[54]=0;i2_ptr[13]=reinterpret_cast<double*>(A);

	// Clear the matrix entries
	for(double *Ap=Am,*Ae=Am+Asize;Ap<Ae;) *(Ap++)=0;

	// Scan the grid and fill in the restriction values
	rt_scan_grid(2048);

	// Send the contributions to other processors
	communicate_rat_strips();

	// If exact solutions are enabled on this level, then perform the LU
	// decomposition of the matrix entries.
	if(Aexact!=NULL) lu_decomposition();
}

/** Performs a scan over the entire grid to do interpolation, restriction, or
 * RAT matrix setup, taking into account the boundary conditions.
 * \param[in] mask an integer that determines the operation to perform. */
void region::rt_scan_grid(unsigned int mask) {
	bool erow=!z_prd&&(p->o&1)==0&&kp==op-1;
	if(p->ak&1) rt_scan_slice(-1,128|512|mask);
	for(int k=0;k<(erow?so-2:so-1);k++) rt_scan_slice(k,64|128|mask);
	if(erow) rt_scan_slice(so-2,64|128|256|512|mask);
	else rt_scan_slice(so-1,(p->bk&1?64|512:64|128|512)|mask);
}

/** Performs a scan over a full xy slice to do interpolation, restriction, or
 * RAT matrix setup, taking into account the boundary conditions.
 * \param[in] k the slice to consider.
 * \param[in] mask an integer that determines the operation to perform. */
void region::rt_scan_slice(int k,unsigned int mask) {
	bool erow=!y_prd&&(p->n&1)==0&&jp==np-1;
	if(p->aj&1) rt_scan_row(-1,k,16|512|mask);
	for(int j=0;j<(erow?sn-2:sn-1);j++) rt_scan_row(j,k,8|16|mask);
	if(erow) rt_scan_row(sn-2,k,8|16|32|512|mask);
	else rt_scan_row(sn-1,k,(p->bj&1?8|512:8|16|512)|mask);
}

/** Performs a scan over a full row in x to do interpolation, restriction, or
 * RAT matrix setup, taking into account the boundary conditions.
 * \param[in] (j,k) the indices of the row to consider.
 * \param[in] mask an integer that determines the operation to perform. */
void region::rt_scan_row(int j,int k,unsigned int mask) {
	bool erow=!x_prd&&(p->m&1)==0&&ip==mp-1;
	if(mask&2048) {
		if(mask&1024) {

			// Perform the RAT bound computation. Depending on
			// whether this row touches any boundary, which is
			// signified by the 512 bit in the mask being set, call
			// the full RAT bound routine or the partial RAT bound
			// routine.
			if(p->ai&1) rat_bound_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_bound_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_bound_full(i,j,k);
			}
			if(erow) rat_bound_partial(sm-2,j,k,(1|2|4)|mask);
			else rat_bound_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		} else {

			// Fill in the RAT matrix entries, again selecting
			// either the full setup routine or the boundary
			// routine depending on the mask value
			if(p->ai&1) rat_fill_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_fill_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_fill_full(i,j,k);
			}
			if(erow) rat_fill_partial(sm-2,j,k,(1|2|4)|mask);
			else rat_fill_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		}
	} else {
		if(mask&1024) {

			// Perform a restriction
			if(p->ai&1) restrict_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) restrict_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) restrict_full(i,j,k);
			}
			if(erow) restrict_partial(sm-2,j,k,(1|2|4)|mask);
			else restrict_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		} else {

			// Perform an interpolation
			if(p->ai&1) interpolate_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) interpolate_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) interpolate_full(i,j,k);
			}
			if(erow) interpolate_partial(sm-2,j,k,(1|2|4)|mask);
			else interpolate_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		}
	}
}

/** Interpolates a 2 by 2 by 2 cube of values in the parent grid, for the case
 * when all of the solution values needed are local to this processor.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to interpolate. */
inline void region::interpolate_full(int i,int j,int k) {

	// Set pointers to source and destination arrays
	double *pp=x0+(i+mg*(j+ng*k)),*pxp=px0+(2*(i+pmg*j+pmng*k)),

	// Construct references to all points in the source array that are
	// involved
	       &v0=*pp,&v1=pp[1],&v2=pp[mg],&v3=pp[mg+1],t1,t2,t3,
	       &v4=pp[mng],&v5=pp[mng+1],&v6=pp[mng+mg],&v7=pp[mng+mg+1];

	// Perform the interpolation to the 2x2x2 block of the destination
	// array, saving several repeated quantities in t1, t2, and t3
	*pxp+=v0;
	pxp[1]+=0.5*(t1=v0+v1);
	pxp[pmg]+=0.5*(v0+v2);
	pxp[pmg+1]+=0.25*(t2=t1+v2+v3);
	pxp[pmng]+=0.5*(v0+v4);
	pxp[pmng+1]+=0.25*(t1+(t3=v4+v5));
	pxp[pmng+pmg]+=0.25*(v0+v2+v4+v6);
	pxp[pmng+pmg+1]+=0.125*(t2+t3+v6+v7);
}

/** Interpolates some part of a 2 by 2 by 2 cube of values in the parent grid,
 * handling the cases when some of the solution values may live on other
 * processors.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to interpolate. */
void region::interpolate_partial(int i,int j,int k,unsigned int mask) {

	// Determine the indices of the points to interpolate from
	int ii=i+1,jj=j+1,kk=k+1;

	// Set pointers the destination array
	double *pxp=px0+2*(i+pmg*j+pmng*k),

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	       &v0=iref(i,j,k),&v1=iref(ii,j,k),&v2=iref(i,jj,k),&v3=iref(ii,jj,k),t2,t3,
	       &v4=iref(i,j,kk),&v5=iref(ii,j,kk),&v6=iref(i,jj,kk),&v7=iref(ii,jj,kk);

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) *pxp+=v0;
	if((mask&74)==74) pxp[1]+=xs?v1:0.5*(v0+v1);
	if((mask&81)==81) pxp[pmg]+=ys?v2:0.5*(v0+v2);
	if((mask&137)==137) pxp[pmng]+=zs?v4:0.5*(v0+v4);
	if((mask&138)==138) pxp[pmng+1]+=zs?(xs?v5:0.5*(v4+v5)):(xs?0.5*(v1+v5):0.25*(v0+v1+v4+v5));
	if((mask&145)==145) pxp[pmng+pmg]+=zs?(ys?v6:0.5*(v4+v6)):(ys?0.5*(v2+v6):0.25*(v0+v2+v4+v6));
	if((mask&18)==18) {
		t2=ys?(xs?v3:0.5*(v2+v3)):(xs?0.5*(v1+v3):0.25*(v0+v1+v2+v3));
		if((mask&82)==82) pxp[pmg+1]+=t2;
		if((mask&146)==146) {
			t3=ys?(xs?v7:0.5*(v6+v7)):(xs?0.5*(v5+v7):0.25*(v4+v5+v6+v7));
			pxp[pmng+pmg+1]+=zs?t3:0.5*(t3+t2);
		}
	}
}

/** Computes the contributions to the restriction from a 2 by 2 by 2 cube of
 * values in the parent grid, for the case when the contributions will be
 * stored local to this processor.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to store the restriction
 *		      contributions. */
inline void region::restrict_full(int i,int j,int k) {

	// Set pointers to source and destination arrays
	double *pp=r0+(i+sm*j+smn*k),&v0=*pp,&v1=pp[1],&v2=pp[sm],&v3=pp[sm+1],
	       &v4=pp[smn],&v5=pp[smn+1],&v6=pp[smn+sm],&v7=pp[smn+sm+1],rd;

	// Compute the grid index in the parent region to consider
	int ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1);

	// Consider the 2 by 2 by 2 set of points in the parent region, compute
	// the residual, and add the contributions to the correct grid points
	// in this region
	v0+=0.125*p->residual(ti,tj,tk);
	rd=0.0625*p->residual(ti+1,tj,tk);v0+=rd;v1+=rd;
	rd=0.0625*p->residual(ti,tj+1,tk);v0+=rd;v2+=rd;
	rd=0.03125*p->residual(ti+1,tj+1,tk);v0+=rd;v1+=rd;v2+=rd;v3+=rd;
	rd=0.0625*p->residual(ti,tj,tk+1);v0+=rd;v4+=rd;
	rd=0.03125*p->residual(ti+1,tj,tk+1);v0+=rd;v1+=rd;v4+=rd;v5+=rd;
	rd=0.03125*p->residual(ti,tj+1,tk+1);v0+=rd;v2+=rd;v4+=rd;v6+=rd;
	rd=0.015625*p->residual(ti+1,tj+1,tk+1);
	v0+=rd;v1+=rd;v2+=rd;v3+=rd;
	v4+=rd;v5+=rd;v6+=rd;v7+=rd;
}

/** Computes the contributions to the restriction from a partial 2 by 2 by 2
 * cube of values in the parent grid, for the case when some contributions might
 * need to be communicated to neighboring processors.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to store the restriction
 *		      contributions. */
void region::restrict_partial(int i,int j,int k,unsigned int mask) {

	// Determine the indices of the points to interpolate from
	int ii=i+1,jj=j+1,kk=k+1,ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1);

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	double &v0=iref(i,j,k),&v1=iref(ii,j,k),&v2=iref(i,jj,k),&v3=iref(ii,jj,k),
	       &v4=iref(i,j,kk),&v5=iref(ii,j,kk),&v6=iref(i,jj,kk),&v7=iref(ii,jj,kk),rd;

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) v0+=0.125*p->residual(ti,tj,tk);
	if((mask&74)==74) {
		rd=0.125*p->residual(ti+1,tj,tk);if(xs) v1+=rd;else {v0+=0.5*rd;v1+=0.5*rd;}
	}
	if((mask&81)==81) {
		rd=0.125*p->residual(ti,tj+1,tk);if(ys) v2+=rd;else {v0+=0.5*rd;v2+=0.5*rd;}
	}
	if((mask&82)==82) {
		rd=0.125*p->residual(ti+1,tj+1,tk);
		if(ys) {
			if(xs) v3+=rd;else {v2+=0.5*rd;v3+=0.5*rd;}
		} else {
			if(xs) {v1+=0.5*rd;v3+=0.5*rd;} else {v0+=0.25*rd;v1+=0.25*rd;v2+=0.25*rd;v3+=0.25*rd;}
		}
	}
	if((mask&137)==137) {
		rd=0.125*p->residual(ti,tj,tk+1);if(zs) v4+=rd;else {v0+=0.5*rd;v4+=0.5*rd;}
	}
	if((mask&138)==138) {
		rd=0.125*p->residual(ti+1,tj,tk+1);
		if(zs) {
			if(xs) v5+=rd;else {v4+=0.5*rd;v5+=0.5*rd;}
		} else {
			if(xs) {v1+=0.5*rd;v5+=0.5*rd;} else {v0+=0.25*rd;v1+=0.25*rd;v4+=0.25*rd;v5+=0.25*rd;}
		}
	}
	if((mask&145)==145) {
		rd=0.125*p->residual(ti,tj+1,tk+1);
		if(zs) {
			if(xs) v6+=rd;else {v4+=0.5*rd;v6+=0.5*rd;}
		} else {
			if(xs) {v2+=0.5*rd;v6+=0.5*rd;} else {v0+=0.25*rd;v2+=0.25*rd;v4+=0.25*rd;v6+=0.25*rd;}
		}
	}
	if((mask&146)==146) {
		rd=0.125*p->residual(ti+1,tj+1,tk+1);
		if(!zs) {
			rd*=0.5;
			if(ys) {
				if(xs) v3+=rd;else {v2+=0.5*rd;v3+=0.5*rd;}
			} else {
				if(xs) {v1+=0.5*rd;v3+=0.5*rd;} else {v0+=0.25*rd;v1+=0.25*rd;v2+=0.25*rd;v3+=0.25*rd;}
			}
		}
		if(ys) {
			if(xs) v7+=rd;else {v6+=0.5*rd;v7+=0.5*rd;}
		} else {
			if(xs) {v5+=0.5*rd;v7+=0.5*rd;} else {v4+=0.25*rd;v5+=0.25*rd;v6+=0.25*rd;v7+=0.25*rd;}
		}
	}
}

/** Computes the box bounds
 * \param[in] (i,j,k) */
void region::rat_bound_full(int i,int j,int k) {

	// Set pointers to source and destination arrays
	int *Sp=S+6*(i+sm*(j+sn*k)),
	    ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak;

	// Deal with first layer
	b_red(ti,tj,tk);b_ex(Sp,ui,uj,uk);
	b_red(ti+1,tj,tk);b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_red(ti,tj+1,tk);b_ex(Sp,ui,uj,uk);b_ex(Sp+6*sm,ui,uj+1,uk);
	b_red(ti+1,tj+1,tk);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_ex(Sp+6*sm,ui,uj+1,uk);b_ex(Sp+6*(sm+1),ui+1,uj+1,uk);

	// Deal with second layer
	b_red(ti,tj,tk+1);b_ex(Sp,ui,uj,uk);b_ex(Sp+6*smn,ui,uj,uk+1);
	b_red(ti+1,tj,tk+1);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_ex(Sp+6*smn,ui,uj,uk+1);b_ex(Sp+6*(smn+1),ui+1,uj,uk+1);
	b_red(ti,tj+1,tk+1);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6*sm,ui,uj+1,uk);
	b_ex(Sp+6*smn,ui,uj,uk+1);b_ex(Sp+6*(sm+smn),ui,uj+1,uk+1);
	b_red(ti+1,tj+1,tk+1);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_ex(Sp+6*sm,ui,uj+1,uk);b_ex(Sp+6*(sm+1),ui+1,uj+1,uk);
	b_ex(Sp+6*smn,ui,uj,uk+1);b_ex(Sp+6*(smn+1),ui+1,uj,uk+1);
	b_ex(Sp+6*(smn+sm),ui,uj+1,uk+1);b_ex(Sp+6*(smn+sm+1),ui+1,uj+1,uk+1);
}

void region::rat_bound_partial(int i,int j,int k,unsigned int mask) {

	// Determine the indices of the points to interpolate from
	// Set pointers to source and destination arrays
	int ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak,ii=i+1,jj=j+1,kk=k+1;

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	int *S0=bref(i,j,k),*S1=bref(ii,j,k),*S2=bref(i,jj,k),*S3=bref(ii,jj,k),
	    *S4=bref(i,j,kk),*S5=bref(ii,j,kk),*S6=bref(i,jj,kk),*S7=bref(ii,jj,kk);

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) {
		b_red(ti,tj,tk);b_ex(S0,ui,uj,uk);
	}
	if((mask&74)==74) {
		b_red(ti+1,tj,tk);b_ex(S1,ui+1,uj,uk);
		if(!xs) b_ex(S0,ui,uj,uk);
	}
	if((mask&81)==81) {
		b_red(ti,tj+1,tk);b_ex(S2,ui,uj+1,uk);
		if(!ys) b_ex(S0,ui,uj,uk);
	}
	if((mask&82)==82) {
		b_red(ti+1,tj+1,tk);b_ex(S3,ui+1,uj+1,uk);
		if(!xs) b_ex(S2,ui,uj+1,uk);
		if(!ys) {b_ex(S1,ui+1,uj,uk);if(!xs) b_ex(S0,ui,uj,uk);}
	}
	if((mask&137)==137) {
		b_red(ti,tj,tk+1);b_ex(S4,ui,uj,uk+1);
		if(!zs) b_ex(S0,ui,uj,uk);
	}
	if((mask&138)==138) {
		b_red(ti+1,tj,tk+1);b_ex(S5,ui+1,uj,uk+1);
		if(!xs) b_ex(S4,ui,uj,uk+1);
		if(!zs) {b_ex(S1,ui+1,uj,uk);if(!xs) b_ex(S0,ui,uj,uk);}
	}
	if((mask&145)==145) {
		b_red(ti,tj+1,tk+1);b_ex(S6,ui,uj+1,uk+1);
		if(!ys) b_ex(S4,ui,uj,uk+1);
		if(!zs) {b_ex(S2,ui,uj+1,uk);if(!ys) b_ex(S0,ui,uj,uk);}
	}
	if((mask&146)==146) {
		b_red(ti+1,tj+1,tk+1);
		if(!zs) {
			b_ex(S3,ui+1,uj+1,uk);
			if(!xs) b_ex(S2,ui,uj+1,uk);
			if(!ys) {b_ex(S1,ui+1,uj,uk);if(!xs) b_ex(S0,ui,uj,uk);}
		}
		b_ex(S7,ui+1,uj+1,uk+1);
		if(!xs) b_ex(S6,ui,uj+1,uk+1);
		if(!ys) {b_ex(S5,ui+1,uj,uk+1);if(!xs) b_ex(S4,ui,uj,uk+1);}
	}
}

/** Calculates the RAT matrix contributions for a 2x2x2 box of points when they
 * are not next to any boundary.
 * \param[in] (i,j,k) the lower corner of the 2x2x2 box to consider. */
void region::rat_fill_full(int i,int j,int k) {

	// Set pointers to source and destination arrays
	double **Ap=A+(i+sm*(j+sn*k));
	int *Sp=S+6*(i+sm*(j+sn*k)),
	    ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak;

	// Deal with first layer
	A_red(ti,tj,tk);A_add(*Ap,Sp,ui,uj,uk);
	A_red(ti+1,tj,tk);A_mul(0.5);A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_red(ti,tj+1,tk);A_mul(0.5);A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);
	A_red(ti+1,tj+1,tk);A_mul(0.25);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);A_add(Ap[sm+1],Sp+6*(sm+1),ui+1,uj+1,uk);

	// Deal with second layer
	A_red(ti,tj,tk+1);A_mul(0.5);A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);
	A_red(ti+1,tj,tk+1);A_mul(0.25);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);A_add(Ap[smn+1],Sp+6*(smn+1),ui+1,uj,uk+1);
	A_red(ti,tj+1,tk+1);A_mul(0.25);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);
	A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);A_add(Ap[sm+smn],Sp+6*(sm+smn),ui,uj+1,uk+1);
	A_red(ti+1,tj+1,tk+1);A_mul(0.125);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);A_add(Ap[sm+1],Sp+6*(sm+1),ui+1,uj+1,uk);
	A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);A_add(Ap[smn+1],Sp+6*(smn+1),ui+1,uj,uk+1);
	A_add(Ap[smn+sm],Sp+6*(smn+sm),ui,uj+1,uk+1);A_add(Ap[smn+sm+1],Sp+6*(smn+sm+1),ui+1,uj+1,uk+1);
}

void region::rat_fill_partial(int i,int j,int k,unsigned int mask) {

	// Determine the indices of the points to interpolate from
	// Set pointers to source and destination arrays
	int ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak,ii=i+1,jj=j+1,kk=k+1;

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	int *S0,*S1,*S2,*S3,*S4,*S5,*S6,*S7;
	double *A0,*A1,*A2,*A3,*A4,*A5,*A6,*A7;
	Abref(A0,S0,i,j,k);Abref(A1,S1,ii,j,k);Abref(A2,S2,i,jj,k);Abref(A3,S3,ii,jj,k);
	Abref(A4,S4,i,j,kk);Abref(A5,S5,ii,j,kk);Abref(A6,S6,i,jj,kk);Abref(A7,S7,ii,jj,kk);

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) {
		A_red(ti,tj,tk);A_add(A0,S0,ui,uj,uk);
	}
	if((mask&74)==74) {
		A_red(ti+1,tj,tk);A_sca(xs);A_add(A1,S1,ui+1,uj,uk);
		if(!xs) A_add(A0,S0,ui,uj,uk);
	}
	if((mask&81)==81) {
		A_red(ti,tj+1,tk);A_sca(ys);A_add(A2,S2,ui,uj+1,uk);
		if(!ys) A_add(A0,S0,ui,uj,uk);
	}
	if((mask&82)==82) {
		A_red(ti+1,tj+1,tk);A_sca(xs,ys);A_add(A3,S3,ui+1,uj+1,uk);
		if(!xs) A_add(A2,S2,ui,uj+1,uk);
		if(!ys) {A_add(A1,S1,ui+1,uj,uk);if(!xs) A_add(A0,S0,ui,uj,uk);}
	}
	if((mask&137)==137) {
		A_red(ti,tj,tk+1);A_sca(zs);A_add(A4,S4,ui,uj,uk+1);
		if(!zs) A_add(A0,S0,ui,uj,uk);
	}
	if((mask&138)==138) {
		A_red(ti+1,tj,tk+1);A_sca(xs,zs);A_add(A5,S5,ui+1,uj,uk+1);
		if(!xs) A_add(A4,S4,ui,uj,uk+1);
		if(!zs) {A_add(A1,S1,ui+1,uj,uk);if(!xs) A_add(A0,S0,ui,uj,uk);}
	}
	if((mask&145)==145) {
		A_red(ti,tj+1,tk+1);A_sca(ys,zs);A_add(A6,S6,ui,uj+1,uk+1);
		if(!ys) A_add(A4,S4,ui,uj,uk+1);
		if(!zs) {A_add(A2,S2,ui,uj+1,uk);if(!ys) A_add(A0,S0,ui,uj,uk);}
	}
	if((mask&146)==146) {
		A_red(ti+1,tj+1,tk+1);
		A_sca(xs,ys,zs);
		if(!zs) {
			A_add(A3,S3,ui+1,uj+1,uk);
			if(!xs) A_add(A2,S2,ui,uj+1,uk);
			if(!ys) {A_add(A1,S1,ui+1,uj,uk);if(!xs) A_add(A0,S0,ui,uj,uk);}
		}
		A_add(A7,S7,ui+1,uj+1,uk+1);
		if(!xs) A_add(A6,S6,ui,uj+1,uk+1);
		if(!ys) {A_add(A5,S5,ui+1,uj,uk+1);if(!xs) A_add(A4,S4,ui,uj,uk+1);}
	}
}

/** Computes the RAT bound for a grid point in the parent region class, and
 * reduces it to the RAT bound for this region, storing it in the temporary
 * bound array.
 * \param[in] (i,j,k) the local indices of the grid point. */
void region::b_red(int i,int j,int k) {
	p->Sbound(i,j,k,St);
	grid_map(*St,St[1],p->m,x_prd);
	grid_map(St[2],St[3],p->n,y_prd);
	grid_map(St[4],St[5],p->o,z_prd);
}

/** Computes the RAT bound for a grid point in the parent region class, and
 * reduces it to the RAT bound for this region.
 * \param[in] (i,j,k) the local indices of the grid point.
 * \param[in] St a pointer to an array of six integers in which to store the
 *		 bound, in global indices. */
void region::A_red(int i,int j,int k) {

	// Compute the bounds of the linear system contributions for this
	// particular gridpoint
	p->Sbound(i,j,k,St);
	St[6]=*St;St[7]=St[1];
	St[8]=St[2];St[9]=St[3];
	St[10]=St[4];St[11]=St[5];
	grid_map(*St,St[1],p->m,x_prd);
	grid_map(St[2],St[3],p->n,y_prd);
	grid_map(St[4],St[5],p->o,z_prd);

	// Ensure the there is enough temporary memory to compute the RAT
	// contribution
	int sz=rat_size(St);
	if(Atmem<sz) {
		delete [] At;
		At=new double[sz];
		Atmem=sz;
	}

	// Loop over all of the box of linear system contributions
	int ii,jj,kk,vi,vj,vk,li,lj,lk,ui,uj,uk,ci,cj,ck,
	    t0=St[7]-St[6],t1=St[9]-St[8],dis=(St[10]*t1+St[8])*t0+St[6];
	double *Ap=At,*Apar=p->A_ref(i,j,k),kfac,jkfac;
	for(vk=St[4];vk<St[5];vk++) {
		ck=igrid_map(lk,uk,St[10],St[11],vk,p->o,z_prd);
		for(vj=St[2];vj<St[3];vj++) {
			cj=igrid_map(lj,uj,St[8],St[9],vj,p->n,y_prd);
			for(vi=*St;vi<St[1];vi++,Ap++) {
				ci=igrid_map(li,ui,St[6],St[7],vi,p->m,x_prd);

				// Assemble this linear system contribution by
				// adding up all relevant entries from the
				// parent linear system
				*Ap=0;
				for(kk=lk;kk<uk;kk++) {
					kfac=kk==ck?1:0.5;
					for(jj=lj;jj<uj;jj++) {
						jkfac=jj==cj?kfac:0.5*kfac;
						for(ii=li;ii<ui;ii++){
							*Ap+=0.125*Apar[(kk*t1+jj)*t0+ii-dis]
							    *(ii==ci?jkfac:0.5*jkfac);
                        }
					}
				}
			}
		}
	}
}

/** Adds contributions to the RAT entries at a particular grid point.
 * \param[in] Ap a pointer to the linear system coefficients for this grid
 *		 point.
 * \param[in] Sp a pointer the
 * \param[in] (i,j,k) the indices to displace the temporary RAT bound by. */
void region::A_add(double *Ap,int *Sp,int i,int j,int k) {
	int vi,vj,vk,t0=Sp[1]-*Sp,t1=Sp[3]-Sp[2];
	double *Atp=At;
	k+=Sp[4];j+=Sp[2];i+=*Sp;
	for(vk=St[4]-k;vk<St[5]-k;vk++) for(vj=St[2]-j;vj<St[3]-j;vj++)
		for(vi=*St-i;vi<St[1]-i;vi++)
			Ap[(vk*t1+vj)*t0+vi]+=*(Atp++);
}

/** Adds contributions to the RAT entries at a particular grid point.
 * \param[in] Ap a pointer to the linear system coefficients for this grid
 *		 point.
 * \param[in] Sp a pointer the */
void region::A_add(double *Av,int *Sv,double *Ap,int *Sp) {
	int vi,vj,vk,t0=Sp[1]-*Sp,t1=Sp[3]-Sp[2];
	double *Avp=Av;
	for(vk=Sv[4]-Sp[4];vk<Sv[5]-Sp[4];vk++) for(vj=Sv[2]-Sp[2];vj<Sv[3]-Sp[2];vj++)
		for(vi=*Sv-*Sp;vi<Sv[1]-*Sp;vi++)
			Ap[(vk*t1+vj)*t0+vi]+=*(Avp++);
}

/** Maps a lower and upper RAT bound in global coordinates from a parent region
 * class into the RAT bounds for this class.
 * \param[in,out] (lo,hi) the lower and upper global RAT bound of the parent
 *			  region class, which are mapped into the bounds for
 *			  this problem upon completion of the routine.
 * \param[in] ma size of the global problem in this coordinate direction.
 * \param[in] prd the periodicity in this coordinate direction. */
void region::grid_map(int &lo,int &hi,int ma,bool prd) {
	if(prd) {
		if(lo<-ma) p_fatal_error("S lower bound too small [prd]",1);
		if(lo>ma-1) p_fatal_error("S lower bound too large [prd]",1);
		if(hi<1) p_fatal_error("S upper bound too small [prd]",1);
		if(hi>2*ma) p_fatal_error("S upper bound too large [prd]",1);
		lo=lo<0?(lo-(ma&1?2:1))/2:lo>>1;
		hi=(hi+((ma&1)&&hi>ma?3:2))>>1;
	} else {
		if(lo<0) p_fatal_error("S lower bound too small [non-prd]",1);
		if(lo>ma-1) p_fatal_error("S lower bound too large [non-prd]",1);
		if(hi<1) p_fatal_error("S upper bound too small [non-prd]",1);
		if(hi>ma) p_fatal_error("S upper bound too large [non-prd]",1);
		lo=(!(ma&1)&&lo==ma-1?lo+1:lo)>>1;
		hi=(hi+2)>>1;
	}
}

/** Calculates the range of indices in the parent class that an index in the
 * current class has a restriction contribution with.
 * \param[out] (lo,hi) the lower and upper global bounds of the parent
 *		       region class.
 * \param[in] (blo,bhi) bounds coming from the extent of the parent stencil,
 *			used to crop the computed bounds to ensure they are in
 *			range.
 * \param[in] i the index in the current class to consider.
 * \param[in] ma size of the global parent problem in this coordinate
 *		 direction.
 * \param[in] prd the periodicity in this coordinate direction.
 * \return The parent index corresponding to the given index. */
int region::igrid_map(int &lo,int &hi,int blo,int bhi,int i,int ma,bool prd) {
	int j;
	if(prd) {
		j=i*2;
		if(ma&1) {
			int llo,hhi;
			if(j<0) {j+=1;llo=-ma;hhi=0;}
			else if(j>=ma) {j-=1;llo=ma;hhi=2*ma;}
			else {llo=0;hhi=ma;}
			lo=j==llo?llo:j-1;
			hi=j==hhi-1?hhi:j+2;
		} else {lo=j-1;hi=j+2;}
	} else {
		if(i<0) p_fatal_error("Grid index too small [non-prd]",1);
		if(i>=(ma+2)>>1) p_fatal_error("Grid index too large [non-prd]",1);
		j=i<<1;
		if(ma&1) {lo=j==0?0:j-1;hi=j==ma-1?ma:j+2;}
		else {
			if(j==ma) {lo=ma-1;hi=ma;j=ma-1;}
			else {
				lo=j==0?0:j-1;
				hi=j==ma-2?ma-1:j+2;
			}
		}
	}
	if(lo<blo) lo=blo;
	if(hi>bhi) hi=bhi;
	return j;
}

/** Assembles the box bound of coefficients for a particular grid point in terms
 * of the global problem indices.
 * \param[in] (i,j,k) the local indices of the grid point to consider.
 * \param[in] Sv a pointer to an array of six integers in which to store the
 *		  bound. */
void region::Sbound(int i,int j,int k,int *Sv) {
	int *Sp=S+6*(i+sm*(j+sn*k));i+=ai;j+=aj;k+=ak;
	*Sv=*Sp+i;Sv[1]=Sp[1]+i;
	Sv[2]=Sp[2]+j;Sv[3]=Sp[3]+j;
	Sv[4]=Sp[4]+k;Sv[5]=Sp[5]+k;
}

/** Extends a RAT bound to contain the current temporary RAT bound.
 * \param[in] Sp the RAT bound to extend.
 * \param[in] (i,j,k) the indices to displace the temporary RAT bound by. */
void region::b_ex(int *Sp,int i,int j,int k) {
	if(*St-i<*Sp) *Sp=*St-i;
    if(St[1]-i>Sp[1]) Sp[1]=St[1]-i;
	if(St[2]-j<Sp[2]) Sp[2]=St[2]-j;
    if(St[3]-j>Sp[3]) Sp[3]=St[3]-j;
	if(St[4]-k<Sp[4]) Sp[4]=St[4]-k;
    if(St[5]-k>Sp[5]) Sp[5]=St[5]-k;
}

/** Extends a RAT bound to contain a given RAT bound.
 * \param[in] Sv the given RAT bound.
 * \param[in] Sp the RAT bound to extend. */
void region::b_ex(int *Sv,int *Sp) {
	if(*Sv<*Sp) *Sp=*Sv;
    if(Sv[1]>Sp[1]) Sp[1]=Sv[1];
	if(Sv[2]<Sp[2]) Sp[2]=Sv[2];
    if(Sv[3]>Sp[3]) Sp[3]=Sv[3];
	if(Sv[4]<Sp[4]) Sp[4]=Sv[4];
    if(Sv[5]>Sp[5]) Sp[5]=Sv[5];
}

/** Returns a reference to the field value at a particular grid index. This may
 * either be within the memory for a primary grid, or within a memory buffer
 * allocated to a ghost region.
 * \param[in] (i,j,k) the grid indices.
 * \return A reference to the field value. */
double& region::iref(int i,int j,int k) {
	int o=(i>=0?(i<sm?1:2):0)+(j>=0?(j<sn?3:6):0)+(k>=0?(k<so?9:18):0),
	    *ixp=i2_inf+o*4;
	return i2_ptr[o]!=NULL?i2_ptr[o][i+*ixp*(j+k*ixp[1])-ixp[2]]:out_of_bounds;
}

/** Returns a pointer to the bound information at a particular grid index.
 * This may either be within the memory for a primary grid, or within a memory
 * buffer allocated to a ghost region.
 * \param[in] (i,j,k) the grid indices.
 * \return A reference to the field value. */
int* region::bref(int i,int j,int k) {
	int o=(i>=0?(i<sm?1:2):0)+(j>=0?(j<sn?3:6):0)+(k>=0?(k<so?9:18):0),
	    *ixp=i2_inf+o*4;
	return S2_ptr[o]!=NULL?S2_ptr[o]+6*(i+*ixp*(j+k*ixp[1])-ixp[2]):NULL;
}

/** Returns pointers to the bound information at a particular grid index.
 * This may either be within the memory for a primary grid, or within a memory
 * buffer allocated to a ghost region.
 * \param[in] (i,j,k) the grid indices.
 * \return A reference to the field value. */
void region::Abref(double *&Ap,int *&Sp,int i,int j,int k) {
	int o=(i>=0?(i<sm?1:2):0)+(j>=0?(j<sn?3:6):0)+(k>=0?(k<so?9:18):0),
	    *ixp=i2_inf+o*4;
	if(i2_ptr[o]==NULL) {
		Ap=NULL;
		Sp=NULL;
	} else {
		Ap=reinterpret_cast<double***>(i2_ptr)[o][i+*ixp*(j+k*ixp[1])-ixp[2]];
		Sp=S2_ptr[o]+6*(i+*ixp*(j+k*ixp[1])-ixp[2]);
	}
}

/** This function fills in the ghost regions of the solution grid with values
 * from neighboring processors. */
void region::communicate() {
	double *pp=com.buf,*outp,**cxp=c_ptr+cneigh_in;
	int i,j,k,*cp=c_inf;
	MPI_Request *reqp=req;

	// Receive the incoming buffers
	for(;cp<c_inf+6*cneigh_in;pp+=cp[5],cp+=6,reqp++)
		MPI_Irecv(pp,cp[5],MPI_DOUBLE,*cp,cp[1],cart,reqp);

	// Prepare and send the outgoing buffers
	for(;cp<c_inf+6*(cneigh_in+cneigh_out);cp+=6,reqp++,cxp++) {
		outp=pp;
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++)
			*(pp++)=(*cxp)[i+mg*(j+ng*k)];
		MPI_Isend(outp,cp[5],MPI_DOUBLE,*cp,cp[1],cart,reqp);
	}

	// Copy the incoming data into the grid
	MPI_Waitall(cneigh_in,req,stat);
	cxp=c_ptr;pp=com.buf;
	for(cp=c_inf;cp<c_inf+6*cneigh_in;cp+=6,cxp++) {
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++)
			(*cxp)[i+mg*(j+ng*k)]=*(pp++);
	}

	// Wait for the outgoing messages to complete
	MPI_Waitall(cneigh_out,req+cneigh_in,stat);
}

/** Communicates any strips of points that are needed during the interpolation
 * step. */
void region::communicate_interpolation_strips() {
	double *pp=com.buf,**ixp=i_ptr,*outp;
	int i,j,k,*cp=i_inf;
	MPI_Request *reqp=req;

	// Receive the incoming buffers
	for(;cp<i_inf+6*ineigh_in;pp+=cp[5],cp+=6,reqp++)
		MPI_Irecv(pp,cp[5],MPI_DOUBLE,*cp,cp[1]|msg_interp,cart,reqp);

	// Prepare and send the outgoing buffers
	for(;cp<i_inf+6*(ineigh_in+ineigh_out);cp+=6,reqp++) {
		outp=pp;
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++)
			*(pp++)=(*ixp)[i+mg*(j+ng*k)];
		MPI_Isend(outp,cp[5],MPI_DOUBLE,*cp,cp[1]|msg_interp,cart,reqp);
		ixp++;
	}
	MPI_Waitall(ineigh_in+ineigh_out,req,stat);
}

/** Communicates any strips of points that are needed during the restriction
 * step. */
void region::communicate_restriction_strips() {
	double *pp=com.buf,**irp=r_ptr,*inp;
	int i,j,k,*cp=i_inf;
	MPI_Request *reqp=req;

	// Send the outgoing buffers
	for(;cp<i_inf+6*ineigh_in;pp+=cp[5],cp+=6,reqp++)
		MPI_Isend(pp,cp[5],MPI_DOUBLE,*cp,cp[1]|msg_rest,cart,reqp);

	// Receive the incoming buffers
	inp=pp;
	for(;cp<i_inf+6*(ineigh_in+ineigh_out);pp+=cp[5],cp+=6,reqp++)
		MPI_Irecv(pp,cp[5],MPI_DOUBLE,*cp,cp[1]|msg_rest,cart,reqp);

	// Add the restriction contributions from the other processors
	MPI_Waitall(ineigh_out,req+ineigh_in,stat);
	for(cp=i_inf+6*ineigh_in;cp<i_inf+6*(ineigh_in+ineigh_out);cp+=6,irp++)
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++)
			(*irp)[i+sm*(j+sn*k)]+=*(inp++);
	MPI_Waitall(ineigh_in,req,stat);
}

/** This communicates the sizes of the RAT matrices in the. */
void region::communicate_rat_bound_strips() {
	int *pp=Strans,*inp,**Spp=S_ptr,i,j,k,*cp=i_inf,l=0,ls=0,*p2=Ac_size;
	MPI_Request *reqp=req;

	// Send the outgoing buffers
	for(;cp<i_inf+6*ineigh_in;cp+=6,reqp++) {
		MPI_Isend(pp,6*cp[5],MPI_INT,*cp,cp[1]|msg_ratb,cart,reqp);
		for(int *pe=pp+6*cp[5];pp<pe;l+=rat_size(pp),pp+=6);
		*(p2++)=l-ls;ls=l;
	}

	// Receive the incoming buffers
	Strans2=inp=pp;
	for(;cp<i_inf+6*(ineigh_in+ineigh_out);pp+=6*cp[5],cp+=6,reqp++)
		MPI_Irecv(pp,6*cp[5],MPI_INT,*cp,cp[1]|msg_ratb,cart,reqp);

	// Add the RAT bound contributions to the relevant gridpoints
	MPI_Waitall(ineigh_out,req+ineigh_in,stat);
	for(cp=i_inf+6*ineigh_in;cp<i_inf+6*(ineigh_in+ineigh_out);cp+=6,Spp++) {
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++,inp+=6) {
			b_ex(inp,(*Spp)+6*(i+sm*(j+sn*k)));
			l+=rat_size(inp);
		}
		*(p2++)=l-ls;ls=l;
	}
	MPI_Waitall(ineigh_in,req,stat);

	// Check that there will be enough space for sending the RAT coefficients
	com.check_buf(l);
}

/** This function fills in the ghost regions of the solution grid with values
 * from neighboring processors. */
void region::communicate_rat_strips() {
	double *pp=com.buf,*inp,***App=A_ptr;
	int **Spp=S_ptr,i,j,k,ijk,*cp=i_inf,*p3=Ac_size,*Sp=Strans2;
	MPI_Request *reqp=req;

	// Send the outgoing buffers
	for(;cp<i_inf+6*ineigh_in;pp+=*p3,cp+=6,p3++,reqp++)
		MPI_Isend(pp,*p3,MPI_DOUBLE,*cp,cp[1]|msg_rat,cart,reqp);

	// Receive the incoming buffers
	inp=pp;
	for(;cp<i_inf+6*(ineigh_in+ineigh_out);pp+=*p3,cp+=6,p3++,reqp++)
		MPI_Irecv(pp,*p3,MPI_DOUBLE,*cp,cp[1]|msg_rat,cart,reqp);

	// Add the RAT bound contributions to the relevant gridpoints
	MPI_Waitall(ineigh_out,req+ineigh_in,stat);
	for(cp=i_inf+6*ineigh_in;cp<i_inf+6*(ineigh_in+ineigh_out);cp+=6,Spp++,App++)
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++,inp+=rat_size(Sp),Sp+=6) {
			ijk=i+sm*(j+sn*k);
			A_add(inp,Sp,(*App)[ijk],(*Spp)+6*ijk);
		}
	MPI_Waitall(ineigh_in,req,stat);
}

/** Computes the global L2 error by evaluating the local L2 error and summing
 * across all other processors at this level.
 * \return The global error if this processor is rank 0, and the local error
 * for all others. */
double region::l2_error_all() {
	int i,j,k;double err=0,rem,sum=0;

	// Update ghost values
	communicate();

	// Compute the local error
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		rem=residual(i,j,k);
		err+=rem*rem;
	}

	// Send all the local results to the zeroth processor
	MPI_Allreduce(&err,&sum,1,MPI_DOUBLE,MPI_SUM,cart);
	return sum;
}

double region::l0_error_all() {
	int i,j,k;double err=0, rem, max=0;

	// Update ghost values
	communicate();

	// Compute the local error
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		rem=fabs(residual(i,j,k));
        if(rem>err) err=rem;
	}

	// Send all the local results to the zeroth processor
	MPI_Allreduce(&err,&max,1,MPI_DOUBLE,MPI_MAX,cart);
	return max;
}


/** Computes the global L2 error by evaluating the local L2 error and summing
 * across all other processors at this level.
 * \return The global error if this processor is rank 0, and the local error
 * for all others. */
double region::l2_error() {
	int i,j,k;double err=0,rem,sum=0;

	// Update ghost values
	communicate();

	// Compute the local error
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		rem=residual(i,j,k);
		err+=rem*rem;
	}

	// Send all the local results to the zeroth processor
	MPI_Reduce(&err,&sum,1,MPI_DOUBLE,MPI_SUM,0,cart);
	return rank==0?sum:err;
}

double region::l0_error() {
	int i,j,k;double err=0, rem, max=0;

	// Update ghost values
	communicate();

	// Compute the local error
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		rem=fabs(residual(i,j,k));
        if(rem>err) err=rem;
	}

	// Send all the local results to the zeroth processor
	MPI_Reduce(&err,&max,1,MPI_DOUBLE,MPI_MAX,0,cart);
	return max;
}

/** Compute the l2 norm, not sqrt, of the rhs */
double region::rhs_l2norm() {
	int i;
	double l_l2=0;
	double g_l2=0;
	for (i =0; i<smno; i++) {
		l_l2+= r0[i]*r0[i];
	}
	MPI_Reduce(&l_l2,&g_l2,1,MPI_DOUBLE,MPI_SUM,0,cart);
	return rank==0?g_l2:l_l2;
}

/** Outputs a z-slice of the x field to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
void region::output_x(const char* filename,int k) {

	// See whether the value of k is within the range of this processor
	int kk=k-ak;
	if(kk>=0&&kk<so) {
		int i,j;float *outp=reinterpret_cast<float*>(com.buf);
		for(j=0;j<sn;j++) for(i=0;i<sm;i++) *(outp++)=x0[i+mg*(j+kk*ng)];
		output_internal(filename,k,true);
	} else output_internal(filename,k,false);
}

/** Outputs a z-slice of the r field to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
void region::output_r(const char* filename,int k) {

	// See whether the value of k is within the range of this processor
	int kk=k-ak;
	if(kk>=0&&kk<so) {
		float *outp=reinterpret_cast<float*>(com.buf);;
		double *rp=r0+kk*smn,*re=rp+smn;
		while(rp!=re) *(outp++)=*(rp++);
		output_internal(filename,k,true);
	} else output_internal(filename,k,false);
}

/** Outputs a z-slice of the residual to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
void region::output_residual(const char* filename,int k) {

	// See whether the value of k is within the range of this processor
	int kk=k-ak;
	if(kk>=0&&kk<so) {
		int i,j;float *outp=reinterpret_cast<float*>(com.buf);
		for(j=0;j<sn;j++) for(i=0;i<sm;i++) *(outp++)=residual(i,j,kk);
		output_internal(filename,k,true);
	} else output_internal(filename,k,false);
}

/** The internal part of the output routine. It assumes that a cross-section of
 * a field has been placed in the out_matrix array. The zeroth processor
 * gathers all of the the local field components and outputs the complete cross
 * section to a file.
 * \param[in] filename the filename to use.
 * \param[in] data whether or not this processor is sending data for printing.
 */
void region::output_internal(const char* filename,int k,bool data) {
	float *outp=reinterpret_cast<float*>(com.buf);

	// Send data if needed
	if(data) MPI_Isend(outp,smn,MPI_FLOAT,0,rank,cart,req);

	// If this is the master processor, then collect data and save it
	if(rank==0) {
		int i,j,jm,jd,kl,o[4];
		float **matp=new float*[mp*np],**matpp=matp,*pp=outp;
		MPI_Request* reqp=req+1;

		// Find which z-layer of processors holds this particular k
		// slice
		o[2]=0;kl=0;
		while(kl+oso[o[2]]<=k) {
			kl+=oso[o[2]++];
			if(o[2]==op) p_fatal_error("Output slice out of range",1);
		}

		// Collect data from all other processors
		for(o[1]=0;o[1]<np;o[1]++) for(*o=0;*o<mp;o[0]++,reqp++) {
			MPI_Cart_rank(cart,o,o+3);
			MPI_Irecv(pp,osm[*o]*osn[o[1]],MPI_FLOAT,o[3],o[3],cart,reqp);
			*(matpp++)=pp;pp+=osm[*o]*osn[o[1]];
		}
		MPI_Waitall(mp*np,req+1,stat+1);

		// Open file and write header line
		FILE *outf=p_safe_fopen(filename,"wb");
		float *fbuf=new float[m+1];
		*fbuf=m;for(i=0;i<m;i++) fbuf[i+1]=i;
		if(fwrite(fbuf,sizeof(float),m+1,outf)!=(size_t) m+1) p_fatal_error("File output error",1);

		// Write field entries line-by-line
		jm=jd=0;
		for(j=0;j<n;j++) {

			// Write header entry
			*fbuf=j;
			if(fwrite(fbuf,sizeof(float),1,outf)!=1) p_fatal_error("File output error",1);

			// Write line
			for(i=0;i<mp;i++) if(fwrite(matp[jd*mp+i]+osm[i]*jm,sizeof(float),osm[i],outf)!=(size_t) osm[i]) p_fatal_error("File output error",1);

			// Update y position markers
			jm++;
			if(jm==osn[jd]) {jd++;jm=0;}
		}

		// Remove temporary memory and close file
		delete [] fbuf;
		delete [] matp;
		fclose(outf);
	}

	// Wait for data to finish sending
	if(data) MPI_Wait(req,stat);
}

/** Clears all of the x field, including the ghost regions
 */
void region::allclear_x_field() {
	int i,j,k;
	for(k=0;k<og;k++) for(j=0;j<ng;j++) for(i=0;i<mg;i++) x[i+mg*(j+k*ng)]=0;
}

/** Clears the central part of the x field, leaving the ghost regions
 * untouched. */
void region::clear_x_field() {
	int i,j,k;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) x0[i+mg*(j+k*ng)]=0;
}

/** Clears the r field. */
void region::clear_r_field() {
	for(double *rp=r0,*re=r0+smno;rp<re;) *(rp++)=0.;
}

/** Copies the r field into the x field. */
void region::copy_r_to_x() {
	int i,j,k;
	double *rp=r0;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) x0[i+mg*(j+k*ng)]=*(rp++);
}

/** Copies the x field into the r field. */
void region::copy_x_to_r() {
	int i,j,k;
	double *rp=r0;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) *(rp++)=x0[i+mg*(j+k*ng)];
}

/** Sets up the array that gives the dimensions of the other processors, needed
 * for output and also for grid transfers. */
void region::setup_output() {
	osm=new int[mp+np+op+3];osn=osm+mp;oso=osn+np;min_sizes=oso+op;
	if(rank==0)  {

		// On rank 0, check that there will be enough space to receive
		// an entire slice of output, and gather the grid dimensions
		// from the neighboring processors.
		com.check_buf_float(m*n);
		gather_sizes();
	} else {

		// On other ranks, check that there will be enough space to
		// send a local slice of output. For processors that are
		// orthogonally aligned with the (0,0,0) processor, send
		// information about the grid dimension.
		com.check_buf_float(smn);
		if(kp==0) {
			if(jp==0) MPI_Send(&sm,1,MPI_INT,0,msg_trans_dims,cart);
			if(ip==0) MPI_Send(&sn,1,MPI_INT,0,msg_trans_dims,cart);
		} else if(ip==0&&jp==0) MPI_Send(&so,1,MPI_INT,0,msg_trans_dims,cart);
	}

	// Broadcast the grid dimensions and minimum sizes to other processors
	MPI_Bcast(osm,mp+np+op+3,MPI_INT,0,cart);
}

/** A routine run by the master processor to gather information about the
 * dimensions of the other regions, needed for output and for setting up grid
 * transfers. The routine also calculates the minimum dimensions in each
 * direction, which can later be used to determine whether a grid transfer
 * should be created. */
void region::gather_sizes() {
	int q[4],&i=*q,&j=q[1],&k=q[2];j=k=0;
	int &msm=*min_sizes,&msn=min_sizes[1],&mso=min_sizes[2];

	// Receive dimensions in the x direction
	msm=*osm=sm;
	for(i=1;i<mp;i++) {
		MPI_Cart_rank(cart,q,q+3);
		MPI_Recv(osm+i,1,MPI_INT,q[3],msg_trans_dims,cart,stat);
		if(osm[i]<msm) msm=osm[i];
	}

	// Receive dimensions in the y direction
	msn=*osn=sn;i=0;
	for(j=1;j<np;j++) {
		MPI_Cart_rank(cart,q,q+3);
		MPI_Recv(osn+j,1,MPI_INT,q[3],msg_trans_dims,cart,stat);
		if(osn[j]<msn) msn=osn[j];
	}

	// Receive dimensions in the z direction
	mso=*oso=so;j=0;
	for(k=1;k<op;k++) {
		MPI_Cart_rank(cart,q,q+3);
		MPI_Recv(oso+k,1,MPI_INT,q[3],msg_trans_dims,cart,stat);
		if(oso[k]<mso) mso=oso[k];
	}
}

/** Prints out contents of the S matrix for diagnostic purposes. */
void region::diagnostic_S() {
	int i,j,k,*Sp;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		Sp=S+6*(i+sm*(j+sn*k));
		printf("%d: (%d,%d,%d) [%d:%d] [%d:%d] [%d:%d]\n",rank,i,j,k,*Sp,Sp[1],Sp[2],Sp[3],Sp[4],Sp[5]);
	}
}

/** Prints out contents of the A matrix for diagnostic purposes. */
void region::diagnostic_A() {
	int i,j,k,l,ss;double *Amp;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
        double tot=0;
		l=i+sm*(j+sn*k);
		ss=rat_size(S+6*l);Amp=A[l];
		printf("%d: (%d,%d,%d) [",rank,i+ai,j+aj,k+ak);
		for(l=0;l<ss-1;l++) {
            printf("%g,",Amp[l]);
            tot+=Amp[l];
        }
#if PGMG_VERBOSE == 3
		int *Sp=S+6*(i+sm*(j+sn*k));
		printf("%g] {%d:%d} {%d:%d} {%d:%d}\n",Amp[l],*Sp,Sp[1],Sp[2],Sp[3],Sp[4],Sp[5]);
#else
		printf("%g] ",Amp[l]);
        tot+=Amp[l];
        printf(" tot %g\n", tot);
#endif
	}
}

void region::diagnostic_x(){
    printf("rhs size (%d %d %d), solution size (%d %d %d)\n", sm, sn, so, mg, ng, og);
	for(int k=0;k<og;k++) for(int j=0;j<ng;j++) for(int i=0;i<mg;i++) {
        int l = i+ j*mg + k*mng;
		printf("(%d, %d, %d): %g\n",i,j,k, x[l]);
		//printf("            : %g\n",r0[l]);
	}

}

/** Sets up an example function on the grid, which can be used to test the
 * code.
 * \param[in] ca the case to consider.
 * \param[in] xfield true if the x field should be set up, false if the r field
 *		     should be set up. */
void region::setup_test(int ca,bool xfield) {
	double *pp;
	for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++) {
		pp=xfield?x0+(i+mg*(j+ng*k)):r0+(i+sm*(j+sn*k));
		switch(ca) {
			case 0: *pp=i+ai;break;
			case 1: *pp=j+aj;break;
			case 2: *pp=k+ak;break;
			case 3: *pp=ak+aj+ai+k+j+i;break;
			case 4: *pp=rank;break;
			case 5: *pp=i==7&&j==7&&k==7?1:0;break;
			case 6: *pp=1;
		}
	}
}

// ########### FUNCTIONS FOR PCG ##############

/** Computes the residual at a particular gridpoint.
 * \param[in] (i,j,k) the gridpoint to consider, in local index
 * \return The residual. */
double region::residual(int i,int j,int k) {

	int di,dj,dk,l=i+sm*(j+sn*k),*Sp=S+6*l;
	double Ax=0,*Amp=A[l];
	for(dk=k+Sp[4];dk<k+Sp[5];dk++) for(dj=j+Sp[2];dj<j+Sp[3];dj++)
		for(di=i+*Sp;di<i+Sp[1];di++){
			Ax+=(*(Amp++))*x0[di+mg*(dj+ng*dk)];
	}
	return r0[l]-Ax;
}

/** Store the residual in some outside array.
 * \param[out] *ext_r ptr to the outside residual array.
 * Note that the output array has dim of x, indented to first non ghost node.
 */
void region::compute_residual(double *ext_r){
	int i, j, k, ind;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		ind = i+j*mg+k*mng;
		ext_r[ind] = residual(i,j,k);
	}
}

/** Apply the matrix A to a particular grid point in an outside array.
 * \param[in] (i,j,k) grid point in local index.
 * \param[in] *in the outside array to use. Must have mnog dims.
 *   *in also must have been indented to first non-ghost node.
 */
double region::mul_A(int i, int j, int k, const double *in){
	int di,dj,dk,l=i+sm*(j+sn*k),*Sp=S+6*l;
	double Ax=0,*Amp=A[l];
	for(dk=k+Sp[4];dk<k+Sp[5];dk++) for(dj=j+Sp[2];dj<j+Sp[3];dj++)
		for(di=i+*Sp;di<i+Sp[1];di++) Ax+=(*(Amp++))*in[di+mg*(dj+ng*dk)];
	return Ax;
}

/** Apply the matrix A to a particular gridpoint.
 * \param[in] (i,j,k) grid point in local index.
 */
double region::mul_A(int i, int j, int k){
	int di,dj,dk,l=i+sm*(j+sn*k),*Sp=S+6*l;
	double Ax=0,*Amp=A[l];
	for(dk=k+Sp[4];dk<k+Sp[5];dk++) for(dj=j+Sp[2];dj<j+Sp[3];dj++)
		for(di=i+*Sp;di<i+Sp[1];di++) Ax+=(*(Amp++))*x0[di+mg*(dj+ng*dk)];
	return Ax;
}

void region::unpadded_Ax(const double * in, double * out){
    // This is only applicable when the entire problem is on the same grid!
	for(int k=0;k<so;k++) {
        for(int j=0;j<sn;j++) {
            for(int i=0;i<sm;i++) {
                int l=i+sm*j+smn*k,*Sp=S+6*l;
                double Ax=0,*Amp=A[l];
            for(int dk=k+Sp[4];dk<k+Sp[5];dk++) {
                for(int dj=j+Sp[2];dj<j+Sp[3];dj++) {
                    for(int di=i+*Sp;di<i+Sp[1];di++) {
                        int ind = non_neg_mod(di, sm) + sm*non_neg_mod(dj, sn) + smn*non_neg_mod(dk, so);
                        Ax+=(*(Amp++))*in[ind];
                    }
                }
            }
                out[l] = Ax;
            }
        }
	}
}

/** Apply the matrix A to the current solution, x.
 * Input and output array must have dim of x,
 * both must have been indented to first non-ghost node.
 * \param[in] * in pointer to vector to apply A.
 * \param[out] * out pointer to vector to store Ad.
 */
void region::Ax(const double *in, double *out){
	int i,j,k;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		out[i+j*mg+k*mng] = mul_A(i,j,k, in);
	}
}

/** Compute the dot product of the source and solution.
 * This is equivalent to the term r(M^-1r) in PCG.
 * We perform dot product in each subdomain, and sum everything up.
 */
double region::source_soln_dot(){
	int ind, indx;
	double prod=0, tot_prod;
	for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
		ind = i+j*sm+k*smn;
		indx = i+j*mg+k*mng;
		prod += x0[indx] * r0[ind];
	}
	MPI_Allreduce(&prod, &tot_prod, 1, MPI_DOUBLE, MPI_SUM, cart);
	return tot_prod;
}

/** Copy the right hand side into the central region of the solution array.
 */
void region::copy_rhs_to_soln(){
	int ind, indx;
	for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
		ind = i+j*sm+k*smn;
		indx = i+j*mg+k*mng;
		x0[indx] = r0[ind];
	}
	communicate();
}

/** If a verbosity level is selected, this prints out messages about the grid
 * allocation.
 * \param[in] suffix a suffix to add to the message, determining the type of grid
 *		     that has been allocated. */
void region::grid_message(const char *suffix) {
#if PGMG_VERBOSE >= 2
	printf("# Rank %d, %d %d %d %d %d %d (m,n,o)=(%d,%d,%d)%s\n",
	       rank,ai,aj,ak,bi,bj,bk,m,n,o,suffix);
#elif PGMG_VERBOSE == 1
	if(rank==0) printf("# Grid <%d,%d,%d> (m,n,o)=(%d,%d,%d)%s\n",
			   mp,np,op,m,n,o,suffix);
#endif
}
