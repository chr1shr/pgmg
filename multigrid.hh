/** \file multigrid.hh
 * \brief Header file for the multigrid class. */

#ifndef PGMG_MULTI_HH
#define PGMG_MULTI_HH

#include "buffer.hh"
#include "common.hh"
#include "geometry.hh"
#include "region.hh"
#include <limits>

/** \brief A class for performing a 3D parallel multigrid solve.
 *
 * This class can solve a 3D parallel multigrid solve of a problem class that
 * is passed to it. The class creates a hierarchy of grids, making use of the
 * region class that represents a subsection of particular grid. This class can
 * carry out multigrid operations by coordinating the region classes to
 * communicate and pass information between each other. */
class multigrid {
	public:

		/**
		 * helper struct for simple implementation of v_cycle solves
		 *
		 * offers predicted number of iterations before convergence
		 * and comparison of provided tolerance to residual
		 *
		 * requires setup:
		 *   - characteristic problem scale for tolerance (required)
		 *   - iteration reduction period (optional)
		 *   - memory parameter for iteration tracking through time (optional)
		 */
		struct mg_guide {

			/** default iteration reduction period */
			static const int DEF_PERIOD=16;
			/** initial moving variance */
			static const int INIT_VAR=1000000;
			/** default memory parameter in units of inverse period */
			static const int DEF_ALPHA=3;
			/** cap on iterations per time step */
			static const int MAX_ITERS=100;
			/** fudge factor in tolerance setting */
			static const int FUDGE_FACTOR=10000;
			/** number of additional iterations to perform after
			 *  tolerance is reached */
			static const int INSURANCE=2;

			/** 
			 * iteration reduction period
			 *
			 * after this many timesteps, the predictor will attempt
			 * to reduce the number of iterations by 1
			 */
			int period;

			/**
			 * iteration index
			 *
			 * this is reduced by at least 1 after each iteration.
			 * the predicted number of iterations is (iter / period)
			 */
			int iter;

			/** whether the moving average is being used for prediction */
			bool use_avg;
			/** moving average and variance of the iteration index */
			double iter_avg,iter_var;
			/** memory parameter alpha ~ 1/P, with P the number of timesteps
			 *  "remembered" by the moving average */
			double alpha;
			/** problem tolerance */
			double tol;
			/** whether the tolerance has been set by a user */
			bool tol_set_manually;

			/** constructs a new mg_guide with defaults */
			mg_guide():use_avg(false),tol_set_manually(false) {
				set_period(DEF_PERIOD);
				alpha = (1.*DEF_ALPHA)/period;
			}

			/**
			 * set tolerance, given characteristic problem scale.
			 * note this corresponds to error at a single node
			 * */
			inline void set_tol(double p_scale) {

				// machine epsilon
				double eps = std::numeric_limits<double>::epsilon();

				// construct tolerance
				tol = p_scale*p_scale*eps*eps*FUDGE_FACTOR;
				tol_set_manually=true;
			}

			/** turn the moving average on or off */
			inline void average_on() {use_avg=true;}
			inline void average_on(double alpha_) {use_avg=true; alpha=alpha_;}
			inline void average_off() {use_avg=false;}

			/**
			 * sets period and resets other values accordingly. Note this
			 * will not adjust the default value of alpha
			 */
			inline void set_period(int period_) {
				period=period_;
				iter_avg=iter=period;
				iter_var=INIT_VAR;
			}

			/** update iteration index given the number of extra
			 *  iterations were required above previous prediction */
			void update(int extras);
			/** return the predicted number of iterations */
			inline int pred_iters() {return iter/period;};
		};

		/** The number of Gauss--Seidel sweeps to perform when going
		 * down the multigrid hierarchy. */
		int v_down;
		/** The number of Gauss--Seidel sweeps to perform on the bottom
		 * level of the multigrid hierarchy. */
		int v_bottom;
		/** The number of Gauss--Seidel sweeps to perform when going
		 * up the multigrid hierarchy. */
		int v_up;
		/** The rank of this processor. */
		int rank;
		/** The number of multigrid levels that have a region on this
		 * processor. */
		int nlevels;
		/** The total number of multigrid levels in the hierarchy. */
		int mlevels;
		/** helper struct with tolerance and interation tracking/prediction */
		mg_guide guide;
        /** Pointer to the geometry class that set this up*/
        geometry *geom;
		/** An array of pointers to regions on this processor. */
		region **reg;
		/** An array of pointers to additional geometries on this
		 * processor. If a region represents a grid transfer, the array
		 * entry will be a pointer to the new grid geometry. Otherwise
		 * the array entry is a null pointer. */
		geometry **geo;
		template<class p_class>
		multigrid(p_class &pr,geometry &gm,comm_buffer &com_);
		~multigrid();
		template<class p_class>
		void setup_matrices(p_class &pr);
		template<class p_class>
		void setup_fields(p_class &pr);
		template<class p_class>
		void setup_r_field(p_class &pr);
		inline void allclear_x_field(){
				reg[0]->allclear_x_field();
		}
		void compute_rats();
		void output_x(const char* filename,int k,int level=0);
		void output_r(const char* filename,int k,int level=0);
		void output_residual(const char* filename,int k,int level=0);
		void setup_hierarchy(geometry &gm);
		void v_cycle(int mode=0);
		void fmg_cycle();
		double l2_error(int level=0);
		double l2_error_all(int level=0);
		double l0_error(int level=0);
		double l0_error_all(int level=0);
		double rhs_l2norm(int level=0);
		void gauss_seidel(int mode=0, int level=0);

		/** 
		 * performs a series of v_cycles until the tolerance in mg_guide
		 * is reached. The number of v_cycles predicted by the guide struct
		 * are performed before the first residual check.
		 *
		 * @param[out] res global l2 error (NOT normalized by number of nodes)
		 * @return 0 if successful convergence, 1 if unsuccessful because
		 *   exceeded maximum iterations defined in mg_pred, 2 if nan present
		 */
		int solve_v_cycles(double &res) {

			// check mg_guide has been set up
			if (!guide.tol_set_manually) p_fatal_error(
				"multigrid::solve_v_cycles cannot"
				" be used before the class mg_guide instance uses"
				" the function multigrid::mg_pred::set_tol to set"
				" the target tolerance to be reached",1);

			// scale tolerance by number of nodes
			double tol = guide.tol * geom->m;
			tol *= geom->n; tol *= geom->o;

			// perform the predicted number of v_cycles
			for (int i=0;i<guide.pred_iters();i++) v_cycle();

			// keep track of extra rounds
			//
			// do 1 more, then 2 more, then 3 more, etc
			int extra_rounds=0,n=0;;
			while ((res=l2_error_all()) > tol && n < guide.MAX_ITERS) {

				// check for nan
				if (std::isnan(res)) return 2;
				
				// do the additional iterations
				extra_rounds++;
				for (int i=0;i<extra_rounds;i++,n++) v_cycle();
			}

			// check if we exceeded cap
			if (n > guide.MAX_ITERS) return 1;

			// do a few extra, update guide, and get out
			for (int i=0;i<guide.INSURANCE;i++) v_cycle();
			guide.update(extra_rounds*(extra_rounds+1)/2);
			return 0;
		}
		int solve_v_cycles() {double res; return solve_v_cycles(res);}

		/** Functions needed for preconditionging conjugate gradient */
		inline void compute_residual(double *ext_r){
			reg[0]->compute_residual(ext_r);
		}
		void Ax(const double *in, double *out){
			reg[0]->Ax(in, out);
		}
		void unpadded_Ax(const double *in, double *out, const int lev=0){
            if(geom->procs!=1 || lev >= nlevels) {
                p_fatal_error("Unpadded_Ax only available for single MPI process, within multigrid hierarchy.", 1);
            }
            else reg[lev]->unpadded_Ax(in, out);
		}
		void copy_rhs_to_soln(){
			return reg[0]->copy_rhs_to_soln();
		}
		double source_soln_dot(){
			return reg[0]->source_soln_dot();
		}
		template<class p_class>
		void update_r_field(p_class &pr){
			reg[0]->update_r_field(pr);
		}
		template<class p_class>
		void overwrite_r_field(p_class &pr){
			reg[0]->overwrite_r_field(pr);
		}
		double soln_diff(double *ext_x){
			return reg[0]->soln_diff(ext_x);
		}

		void interpolation(int level);
		void restriction(int level);
		void communicate(int level=0);
		inline void diagnostic_S(int level=0) {
			if(level<nlevels) reg[level]->diagnostic_S();
		}
		inline void diagnostic_A(int level=0) {
			if(level<nlevels) reg[level]->diagnostic_A();
		}
		inline void diagnostic_x(int level=0) {
			if(level<nlevels) reg[level]->diagnostic_x();
		}
		inline void solve_exact(int level=0, bool smooth=0) {
			if(level<nlevels) reg[level]->solve_exact(smooth);
		}
		inline void free() {
			for(int i=nlevels-1;i>=0;i--) if(geo[i]!=NULL) geo[i]->free();
		}
		inline void setup_test(int ca,bool xfield,int level=0) {
			if(level<nlevels) reg[level]->setup_test(ca,xfield);
		}
	protected:
		comm_buffer &com;
};

/** Initializes the multigrid class, setting up a hierarchy of region classes to handle
 * computations on progressively coarser grids.
 * \param[in] pr a pointer to the problem class to solve.
 * \param[in] comm_ a reference to a communication buffer class. */
template<class p_class>
multigrid::multigrid(p_class &pr,geometry &gm,comm_buffer &com_)
	: v_down(2), v_bottom(20), v_up(2), rank(gm.rank),
	nlevels(0), geom(&gm), com(com_) {

	// Allocate the top level of the hierarchy using the problem-specific
	// setup routines
	reg=new region*[mg_max_levels];
	*reg=new region(pr,gm,com);

	// Set up the rest of the hierarchy. This part is common to all
	// problems and does not need to be templated.
	setup_hierarchy(gm);
}

/** Sets up the matrices on each level, by filling in those on the top level
 * from a problem class, and then computing the RAT matrices for the lower
 * levels.
 * \param[in] pr a pointer to a problem class to use. */
template<class p_class>
void multigrid::setup_matrices(p_class &pr) {
	reg[0]->setup_matrices(pr);
	compute_rats();
}

/** Sets up the source field on the top level from a problem class.
 * \param[in] pr a pointer to a problem class to use. */
template<class p_class>
void multigrid::setup_r_field(p_class &pr) {
	reg[0]->setup_r_field(pr);
}

/** Sets up the fields on the top level from a problem class.
 * \param[in] pr a pointer to a problem class to use. */
template<class p_class>
void multigrid::setup_fields(p_class &pr) {
	reg[0]->setup_fields(pr);
}

#endif
