#ifndef PGMG_GEOMETRY_HH
#define PGMG_GEOMETRY_HH

#include "common.hh"

class region;

class geometry {
	public:
		/** The global x grid size. */
		const int m;
		/** The global y grid size. */
		const int n;
		/** The global z grid size. */
		const int o;
		/** The global x periodicity. */
		const bool x_prd;
		/** The global y periodicity. */
		const bool y_prd;
		/** The global z periodicity. */
		const bool z_prd;
		/** The rank of the processor. */
		int rank;
		/** The total number of processors. */
		int procs;
		/** The number of processors in the x direction. */
		int mp;
		/** The number of processors in the y direction. */
		int np;
		/** The number of processors in the z direction. */
		int op;
		/** The x index of this processor. */
		int ip;
		/** The y index of this processor. */
		int jp;
		/** The z index of this processor. */
		int kp;
		/** The global lower x index (inclusive) for this processor. */
		int ai;
		/** The global lower y index (inclusive) for this processor. */
		int aj;
		/** The global lower z index (inclusive) for this processor. */
		int ak;
		/** The global upper x index (exclusive) for this processor. */
		int bi;
		/** The global upper y index (exclusive) for this processor. */
		int bj;
		/** The global upper z index (exclusive) for this processor. */
		int bk;
		/** The local x grid size. */
		int sm;
		/** The local y grid size. */
		int sn;
		/** The local z grid size. */
		int so;
		/** The product of the local x and y grid sizes, used for
		 * stepping through memory. */
		int smn;
		/** The total number of local grid points. */
		int smno;
		/** The cartesian communicator for this grid. */
		MPI_Comm cart;
		/** The table of processor ranks in the global problem, if
		 * needed. */
		int *ranks;
		geometry(int m_,int n_,int o_,bool x_prd_,bool y_prd_,bool z_prd_);
		geometry(region *reg,int pprocs,int thin);
		geometry(const geometry &gm);
		~geometry();
		void free();
		void set_up_neighbors(int *neigh);
	private:
		void search_optimal(long &msize);
		void create_cart_communicator(MPI_Comm comm);
};

#endif
