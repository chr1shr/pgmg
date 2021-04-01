/** \file region.hh
 * \brief Header file for the region class, representing part of a particular
 * grid level within a multigrid hierarchy. */

#ifndef PGMG_REGION_HH
#define PGMG_REGION_HH

#include "buffer.hh"
#include "common.hh"
#include "geometry.hh"

class region {
	public:
		/** A pointer to a parent region class. */
		region *p;
        /** A static counter for GS sweeps */
        static int gs_num;
		/** The rank of the processor that created this class. */
		int rank;
		/** The periodicity in the x direction. */
		const bool x_prd;
		/** The periodicity in the y direction. */
		const bool y_prd;
		/** The periodicity in the z direction. */
		const bool z_prd;
		/** A flag to enanble Gauss--Seidel sweeps on this level, used
		 * to skip redundant grids in a multigrid hierarchy. */
		bool gs_enable;
		/** The global size of the problem in the x direction. */
		const int m;
		/** The global size of the problem in the y direction. */
		const int n;
		/** The global size of the problem in the z direction. */
		const int o;
		/** The x index of this region within the grid. */
		const int ip;
		/** The y index of this region within the grid. */
		const int jp;
		/** The z index of this region within the grid. */
		const int kp;
		/** The size of the region grid in the first coordinate. */
		const int mp;
		/** The size of the region grid in the second coordinate. */
		const int np;
		/** The size of the region grid in the third coordinate. */
		const int op;
		/** A pointer to the array holding the solution vector,
		 * including ghost regions that are filled with values from
		 * other processors. */
		double *x;
		/** A pointer to the (0,0,0) solution entry in the local
		 * indexing system. */
		double *x0;
		/** A pointer to the source term array. */
		double *r0;
		/** A pointer to the matrix entries. */
		double **A;
		/** A pointer to the matrix entries. */
		double *Am;
		/** The total number of matrix entries. */
		int Asize;
		/** The size of the solution array in the x direction,
		 * including ghost points. */
		int mg;
		/** The size of the solution array in the y direction,
		 * including ghost points. */
		int ng;
		/** The size of the solution array in the z direction,
		 * including ghost points. */
		int og;
		/** The size of an xy slice in the solution array, including
		 * ghost points. */
		int mng;
		/** THe size of indent to take solution to the first non-ghost node.
		*/
		int indent;
		/** An array of neighboring processor IDs in a 3 by 3 by 3 cube
		 * surrounding this processor. */
		int neighbor[27];
		/** The minimum number of gridpoints owned by any processor, in
		 * each of three directions. */
		int *min_sizes;
		/** The MPI communicator, containing all processors that are
		 * involved in this multigrid level. */
		MPI_Comm cart;
		MPI_Comm pcart;
		template<class p_class>
		region(p_class &pr,geometry &gm,comm_buffer &com_);
		region(region* reg,geometry &gm,comm_buffer &com_);
		region(region *p_,comm_buffer &com_);
		~region();
		/** Sets up the matrices for a region class on the top level,
		 * by filling them in from the problem class.
		 * \param[in] pr a pointer to a problem class to use. */
		template<class p_class>
		void setup_matrices(p_class &pr) {
			int i,j,k;
			double *Amp=Am;
			for(k=ak;k<bk;k++) for(j=aj;j<bj;j++)
				for(i=ai;i<bi;i++) if(pr.internal(i,j,k)) pr.fill_entries(i,j,k,Amp);
		}
		/** Sets up the x and r fields for a region class on the top
		 * level, by filling them in from the problem class.
		 * \param[in] pr a pointer to a problem class to use. */
		template<class p_class>
		void setup_fields(p_class &pr) {
			setup_x_field(pr);
			setup_r_field(pr);
		}
		template<class p_class>
		void setup_x_field(p_class &pr) {
			int i,j,k;

			// Set up the x field
			for(k=0;k<so;k++) for(j=0;j<sn;j++)
				for(i=0;i<sm;i++) x0[i+mg*(j+k*ng)]=pr.x_field(i+ai,j+aj,k+ak);
		}
		template<class p_class>
		void setup_r_field(p_class &pr) {
			int i,j,k;

			// Set up the r field
			double *rp=r0;
			for(k=0;k<so;k++) for(j=0;j<sn;j++)
				for(i=0;i<sm;i++) *(rp++)=pr.r_field(i+ai,j+aj,k+ak);
		}
		void clear_r_field();
		void clear_x_field();
		void allclear_x_field();
		void copy_r_to_x();
		void copy_x_to_r();
		void gauss_seidel(int mode=0);
		void enable_exact();
		void solve_exact(bool smooth=0);
		void restriction();
		void interpolation();
		void compute_rats();
		void communicate();
		void fill_rat_bounds();
		double l2_error();
		double l2_error_all();
		double l0_error();
		double l0_error_all();
		double rhs_l2norm();
		void setup_outgoing_table(geometry &gm);
		void setup_incoming_table();
		void setup_incoming_pointers();
		void allocate_transfer_memory(int tr_size_);
		void send_r();
		void send_x();
		void send_S();
		void send_A();
		/** Waits for a grid transfer sending operation to complete. */
		inline void send_wait() {
			MPI_Waitall(tr_size,req,stat);
		}
		void receive_r();
		void receive_x();
		void receive_S();
		void receive_A();
		/** For the case when this region is the outgoing part of a
		 * grid transfer operation, but there is incoming grid transfer
		 * operation on this processor, this routine checks that the
		 * communication buffers will be large enough for all of the
		 * transfer messages. */
		inline void check_transfer_buf() {
			com.check_buf(Asize);
			com.check_buf_int(6*smno);
			com.check_buf(smno);
		}
		void output_x(const char* filename,int k);
		void output_r(const char* filename,int k);
		void output_residual(const char* filename,int k);
		void output_internal(const char* filename,int k,bool data);
		void setup_test(int ca,bool xfield);
		void diagnostic_S();
		void diagnostic_A();
        void diagnostic_x();
		void diagnostic_sums();
		inline int min_size() {
			return min(min_sizes[2],min(min_sizes[1],*min_sizes));
		}
		inline int total() {
			return m*n*o;
		}
		inline void quick_print() {
			int i,j;
			for(j=0;j<ng;j++) for(i=0;i<mg;i++) {
				printf("DI %d %d %g\n",i,j,x[i+mg*(j+10*ng)]);
			}
		}

		/** Functions needed for preconditioning conjugate graident */
		void Ax(const double *in, double *out);
		void unpadded_Ax(const double *in, double *out);
        int non_neg_mod(const int a, const int div){
            int rem = a%div;
            if(rem<0) rem+=div;
            return rem;
        }
		void compute_residual(double *ext_r);
		double residual(int i,int j,int k);
		double mul_A(int i,int j,int k);
		double mul_A(int i,int j,int k, const double *in);
		double source_soln_dot();
		void copy_rhs_to_soln();
		template<class p_class>
		void overwrite_r_field(p_class &pr) {
			int i,j,k;

			// Set up the r field
			double *rp=r0;
			for(k=0;k<so;k++) for(j=0;j<sn;j++)
				for(i=0;i<sm;i++) *(rp++)=pr.reset_r_field(i+ai,j+aj,k+ak);
		}
		template<class p_class>
		void update_r_field(p_class &pr) {
			int i,j,k;

			// Set up the r field
			double *rp=r0;
			for(k=0;k<so;k++) for(j=0;j<sn;j++)
				for(i=0;i<sm;i++) *(rp++)+=pr.update_r_field(i+ai,j+aj,k+ak);
		}
		// Compute the difference btw two solutions
		double soln_diff(double *ext_x){
			int i, j, k, ind;
			double diff=0, tot_diff;
			for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
				ind = i+j*mg+k*mng;
				double tmp=(ext_x[ind] - x0[ind])*(ext_x[ind] - x0[ind]);
				diff+=tmp;
			}
			MPI_Allreduce(&diff, &tot_diff, 1, MPI_DOUBLE, MPI_SUM, cart);
			return tot_diff;
		}
	//protected:
		/** The total number of neighboring processors. */
		int tneigh;
		/** The total number of neighboring processors to send ghost
		 * solution vector points to. */
		int cneigh_in;
		/** The total number of neighboring processors to receive ghost
		 * solution vector points from. */
		int cneigh_out;
		/** The total number of neighboring processors to send
		 * interpolation ghost points to. */
		int ineigh_in;
		/** The total number of neighboring processors to receive
		 * interpolation ghost points from. */
		int ineigh_out;
		/** A pointer to the box bound information. Each grid point has
		 * six integers giving the lower and upper bounds on terms in
		 * the linear system in each of the three coordinate
		 * directions. */
		int *S;
		/** Temporary space for assemble RAT matrix elements. */
		double *At;
		/** The size of the temporary space for RAT matrix elements. */
		int Atmem;
		/** Temporary space for the RAT bounds. */
		int St[12];
		/** The global lower x bound for this processor. */
		int ai;
		/** The global lower y bound for this processor. */
		int aj;
		/** The global lower z bound for this processor. */
		int ak;
		/** The global upper x bound for this processor. */
		int bi;
		/** The global upper y bound for this processor. */
		int bj;
		/** The global upper z bound for this processor. */
		int bk;
		/** The size of the grid on this processor in the x direction.
		 */
		int sm;
		/** The size of the grid on this processor in the y direction.
		 */
		int sn;
		/** The size of the grid on this processor in the z direction.
		 */
		int so;
		/** The size of an xy grid slice on this processor. */
		int smn;
		/** The total number of real gridpoints on this processor. */
		int smno;
		/** The global lower x bound for this processor, including
		 * ghost regions. */
		int li;
		/** The global lower y bound for this processor, including
		 * ghost regions. */
		int lj;
		/** The global lower z bound for this processor, including
		 * ghost regions. */
		int lk;
		/** The global upper x bound for this processor, including
		 * ghost regions. */
		int hi;
		/** The global upper y bound for this processor, including
		 * ghost regions. */
		int hj;
		/** The global upper z bound for this processor, including
		 * ghost regions. */
		int hk;
		int lip,ljp,lkp;
		int hip,hjp,hkp;
		/** The size of the parent region class's solution array in the
		 * x direction, including ghost points. */
		int pmg;
		/** The size of the parent region class's solution array in the
		 * x direction, including ghost points. */
		int png;
		/** The size of the parent region class's xy slice, including
		 * ghost points. */
		int pmng;
		/** The size of the strip transfer buffer. */
		int ntrans;
		/** A pointer to the (0,0,0) solution gridpoint in the parent
                * region class. */
		double *px0;
		/** A pointer to the (0,0,0) RAT bound gridpoint in the parent
		 * region class. */
		int *pS0;
		/** An array of x dimension sizes for the processors in the
		 * global processor configuration, used for output. This is
		 * only initialized on the master processor. */
		int *osm;
		/** An array of y dimension sizes for the processors in the
		 * global processor configuration, used for output. This is
		 * only initialized on the master processor. */
		int *osn;
		/** An array of z dimension sizes for the processors in the
		 * global processor configuration, used for output. This is
		 * only initialized on the master processor. */
		int *oso;
		/** Information about the main communication buffers, used
		 * to fill the ghost regions required for the Gauss-Seidel
		 * smoothing steps. This consists of records of six integers
		 * each, containing the processor ID to communicate with, and
		 * the dimensions of the buffer. */
		int *c_inf;
		/** An array of pointers to the lower corners of the main
		 * communication buffers. */
		double **c_ptr;
		/* Information about the strip communication buffers, used in
		 * the restriction, interpolation, and RAT computation steps.
		 * This consists of records of six integers each, containing
		 * the processor ID to communicate with, and the dimensions of
		 * the buffer. */
		int *i_inf;
		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the x array, used in the
		 * interpolation step. */
		double **i_ptr;
		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the r array, used in the
		 * restriction step. */
		double **r_ptr;
		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the S array, used in the RAT
		 * bound computation step. */
		int **S_ptr;
		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the A array, used in the RAT
		 * computation step. */
		double ***A_ptr;
		/** Sizes of the RAT coefficient strip transfers. */
		int *Ac_size;
		/** Tables used during the restriction, interpolation, and RAT
		 * computation steps to map references to external grid points
		 * to communication buffers. */
		int *i2_inf;
		/** General purpose pointers that are used in the restriction,
		 * interpolation, and RAT computation steps to mark positions
		 * in the communication buffer corresponding to each
		 * received/sent message. */
		double **i2_ptr;
		/** Pointers to the linear system bound memory in the transfer
		 * strips. */
		int **S2_ptr;
		/** A pointer to memory for storing the RAT bound information
		 * that is passed in the strip communication. */
		int *Strans;
		int *Strans2;
		/** Pointers to entries in the communication buffer for sending
		 * RAT contributions. */
		double **Atrans;
		int tr_size;
		int tr_pAsize;
		int tr_psmno;
		int *tr_inf;
		double **tr_x;
		double **tr_r;
		double ***tr_A;
		int **tr_S;
		int *ipiv;
		double *Aexact;
		/** An array of MPI requests. */
		MPI_Request *req;
		/** An array of MPI statuses. */
		MPI_Status *stat;
		/** The size of the MPI request and status arrays. */
		int req_stat_size;
		/** A reference to a common buffer to be used as temporary
		 * space during communications. */
		comm_buffer &com;
		/** A dummy value to be passed as a reference when looking up . */
		double out_of_bounds;
		int transfer_buffer_sizes(int *ctmp,int xd_size,int xu_size,int yd_size,int yu_size,int zd_size,int zu_size,int &n_in,int &n_out);
		void setup_communication();
		void setup_gs_buffers();
		void setup_strip_buffers();
		void setup_rt_pointers();
		void communicate_interpolation_strips();
		void communicate_restriction_strips();
		void communicate_rat_bound_strips();
		void communicate_rat_strips();
		void rt_scan_grid(unsigned int mask);
		void rt_scan_slice(int k,unsigned int mask);
		void rt_scan_row(int j,int k,unsigned int mask);
		inline void interpolate_full(int i,int j,int k);
		void interpolate_partial(int i,int j,int k,unsigned int mask);
		inline void restrict_full(int i,int j,int k);
		void restrict_partial(int i,int j,int k,unsigned int mask);
		void rat_bound_full(int i,int j,int k);
		void rat_bound_partial(int i,int j,int k,unsigned int mask);
		void rat_fill_full(int i,int j,int k);
		void rat_fill_partial(int i,int j,int k,unsigned int mask);
		inline int rat_size(int *Sp) {
			return (Sp[1]-*Sp)*(Sp[3]-Sp[2])*(Sp[5]-Sp[4]);
		}
		double &iref(int i,int j,int k);
		int* bref(int i,int j,int k);
		void Abref(double *&Ap,int *&Sp,int i,int j,int k);
		void gather_sizes();
		void setup_output();
		inline int min(int a,int b) {return a>b?b:a;}
		inline int max(int a,int b) {return a<b?b:a;}
		inline void six_copy(int *&cp2,int *&cp) {
			*(cp2++)=*(cp++);*(cp2++)=*(cp++);
			*(cp2++)=*(cp++);*(cp2++)=*(cp++);
			*(cp2++)=*(cp++);*(cp2++)=*(cp++);
		}
		void b_ex(int *Sp,int i,int j,int k);
		void b_ex(int *Sv,int *Sp);
		void Sbound(int i,int j,int k,int *Sv);
		inline double* A_ref(int i,int j,int k) {
			return A[i+sm*(j+sn*k)];
		}
		void b_red(int i,int j,int k);
		void A_red(int i,int j,int k);
		void A_add(double *Ap,int *Sp,int i,int j,int k);
		void A_add(double *Av,int *Sv,double *Ap,int *Sp);
		void grid_map(int &lo,int &hi,int ma,bool prd);
		int igrid_map(int &lo,int &hi,int blo,int bhi,int i,int ma,bool prd);
		void incoming_index_range(int *os,int *ws,int vp,int a,int b,int &l,int &u);
		void lu_decomposition();
		void grid_message(const char* suffix);
		inline void A_sca(bool b1,bool b2=true,bool b3=true) {
			if(!(b1&&b2&&b3)) A_mul((b1?1:0.5)*(b2?1:0.5)*(b3?1:0.5));
		}
		inline void A_mul(double fac) {
			for(double *Ap=At,*Ae=At+rat_size(St);Ap<Ae;Ap++) *Ap*=fac;
		}
		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
};

/** This constructor sets up a region for a Gauss-Seidel computation, or a
 * top-level multigrid computation.
 * \param[in] *pr a pointer to the problem class to solve.
 * \param[in] rank_ the current processor's rank.
 * \param[in] com_ a reference to a communication buffer class. */
template<class p_class>
region::region(p_class &pr,geometry &gm,comm_buffer &com_)
	: p(NULL), rank(gm.rank), x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
	gs_enable(true), m(gm.m), n(gm.n), o(gm.o),
	ip(gm.ip), jp(gm.jp), kp(gm.kp), mp(gm.mp), np(gm.np), op(gm.op),
	cart(gm.cart), pcart(MPI_COMM_NULL), At(NULL), ai(gm.ai), aj(gm.aj), ak(gm.ak),
	bi(ai+gm.sm), bj(aj+gm.sn), bk(ak+gm.so),
	sm(gm.sm), sn(gm.sn), so(gm.so), smn(gm.smn), smno(gm.smno), c_inf(NULL),
	c_ptr(NULL), i_inf(NULL), i_ptr(NULL), tr_size(0), tr_inf(NULL), Aexact(NULL),
	com(com_) {

	grid_message(" [top]");

	// Set up neighbor table
	gm.set_up_neighbors(neighbor);
	setup_communication();

	// Set up problem storage space
	int i,j,k,r,*Sp=(S=new int[6*smno]);
	A=new double*[smno];
	Asize=0;

	// Compute amount of memory for matrices and ghost grid size
	li=ai;hi=bi;lj=aj;hj=bj;lk=ak;hk=bk;
	for(k=ak;k<bk;k++) for(j=aj;j<bj;j++) for(i=ai;i<bi;i++) {
		r=pr.range_xd(i,j,k);if(i+r<li) li=i+r;*(Sp++)=r;
		r=pr.range_xu(i,j,k);if(i+r>hi) hi=i+r;*(Sp++)=r;
		r=pr.range_yd(i,j,k);if(j+r<lj) lj=j+r;*(Sp++)=r;
		r=pr.range_yu(i,j,k);if(j+r>hj) hj=j+r;*(Sp++)=r;
		r=pr.range_zd(i,j,k);if(k+r<lk) lk=k+r;*(Sp++)=r;
		r=pr.range_zu(i,j,k);if(k+r>hk) hk=k+r;*(Sp++)=r;
		if(pr.internal(i,j,k)) Asize+=pr.mem_size(i,j,k);
	}
	mg=hi-li;ng=hj-lj;og=hk-lk;mng=mg*ng;

	// Allocate memory for problem entries, and set up pointers
	Am=new double[Asize];
	double *Amp=Am,**Ap=A;
	for(k=ak;k<bk;k++) for(j=aj;j<bj;j++) for(i=ai;i<bi;i++,Ap++) {
		if(pr.internal(i,j,k)) {
			*Ap=Amp;
			Amp+=pr.mem_size(i,j,k);
		} else {
			*Ap=pr.external_ptr(i,j,k);
		}
	}

	// Allocate function and source arrays and set up communication buffers
	x=new double[mng*og];
	r0=new double[smno];
	setup_gs_buffers();

	// Check size for output matrix and l2_error
//	osm=new int[32];
	setup_output();
}

#endif
