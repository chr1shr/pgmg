/** \file test_problem.hh
 * \brief Header file for the test_problem class. */

#ifndef PGMG_TEST_PROBLEM_HH
#define PGMG_TEST_PROBLEM_HH

#include "geometry.hh"

/** This class encapsulates all of the routines that can be used to set up the
 * simple problem, including grid size, and the stencils at each grid point. */
class test_problem {
	public:
		/** The periodicity in the x direction. */
		const bool x_prd;
		/** The periodicity in the y direction. */
		const bool y_prd;
		/** The periodicity in the z direction. */
		const bool z_prd;
		/** The global size of the problem in the x direction. */
		const int m;
		/** The global size of the problem in the y direction. */
		const int n;
		/** The global size of the problem in the z direction. */
		const int o;
		int ai,aj,ak,bi,bj,bk,sm,sn,so;
		/** The value of the diagonal term in an identity row,
		 * corresponding to a gridpoint that is held fixed according to
		 * a Dirichlet boundary condition. */
		static double fix[1];
		static double stencil[27];
		/* unknown vector */
		test_problem(geometry &gm) : x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
			m(gm.m), n(gm.n), o(gm.o),
			ai((m*gm.ip)/gm.mp), aj((n*gm.jp)/gm.np), ak((o*gm.kp)/gm.op),
			bi((m*(gm.ip+1))/gm.mp), bj((n*(gm.jp+1))/gm.np), bk((o*(gm.kp+1))/gm.op),
			sm(bi-ai), sn(bj-aj), so(bk-ak) {
		}
		inline int range_xd(int i,int j,int k) {return interior(i,j,k)?-1:0;}
		inline int range_xu(int i,int j,int k) {return interior(i,j,k)?2:1;}
		inline int range_yd(int i,int j,int k) {return interior(i,j,k)?-1:0;}
		inline int range_yu(int i,int j,int k) {return interior(i,j,k)?2:1;}
		inline int range_zd(int i,int j,int k) {return interior(i,j,k)?-1:0;}
		inline int range_zu(int i,int j,int k) {return interior(i,j,k)?2:1;}
		inline int mem_size(int i,int j,int k) {
			p_fatal_error("No need to measure memory with this problem type",1);
			return 0;
		}
		void fill_entries(int i,int j,int k,double *&en);
		inline double x_field(int i,int j,int k) {return 0;}
		inline double r_field(int i,int j,int k) {
			if (i==0 || j==0 || k==0 || i==m-1 || j==n-1 || k==o-1) return 42*fix[0];

			if (i>20&&j>20&&k>20&&i<40&&j<40&&k<40) return 1;

			return 0;
		}
		inline bool internal(int i,int j,int k) {return false;}
		inline double* external_ptr(int i,int j,int k) {
			return interior(i,j,k)?stencil:fix;}
		inline bool interior(int i,int j,int k) {
			return (x_prd||(i!=0&&i!=m-1))&&
			       (y_prd||(j!=0&&j!=n-1))&&
			       (z_prd||(k!=0&&k!=o-1));
		}
};

#endif
