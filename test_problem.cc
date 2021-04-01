/** \file test_problem.cc
 * \brief Function implementations for the test_problem class. */

#include "common.hh"
#include "test_problem.hh"

/** Fills the entries in the table with stencil values for a particular gridpoint.
 * \param[in] (i,j,k) the index of the gridpoint to consider.
 * \param[in,out] en a reference to the pointer to the stencil table. */
void test_problem::fill_entries(int i,int j,int k,double *&en) {
	p_fatal_error("No need to fill entries with this problem type",1);
}

/** An array containing the basic finite-difference stencil for the Laplacian
 * operator. */
//double test_problem::stencil[27]={-1,-1,-1,-1,0,-1,-1,-1,-1,-1,0,-1,0,6,0,-1,0,-1,-1,-1,-1,-1,0,-1,-1,-1,-1};
// Nodal Poisson problem
// double test_problem::stencil[27]={-1,-2,-1,-2,0,-2,-1,-2,-1,-2,0,-2,0,32,0,-2,0,-2,-1,-2,-1,-2,0,-2,-1,-2,-1};
//double test_problem::stencil[27]={1/12.,2/12.,1/12.,2/12.,0,2/12.,1/12.,2/12.,1/12.,2/12.,0,2/12.,0,-32/12.,0,2/12.,0,2/12.,1/12.,2/12.,1/12.,2/12.,0,2/12.,1/12.,2/12.,1/12.};
// Cell center Poisson problem
double test_problem::stencil[27]={0,0,0,0,1,0,0,0,0,0,1,0,1,-6,1,0,1,0,0,0,0,0,1,0,0,0,0};
//double test_problem::stencil[27]={0,0,0,0,0,0,0,0,0,
//					0,0,0,1,-2,1,0,0,0,
//					0,0,0,0,0,0,0,0,0};
//const double test_problem::stencil[27]={1,-2,1,-2,4,-2,1,-2,1,
//					-2,4,-2,4,-8,4,-2,4,-2,
//					1,-2,1,-2,4,-2,1,-2,1};
//const double test_problem::stencil[27]={-1,-1,-1,-1,-1,-1,-1,-1,-1,
//					 -1,-1,-1,-1,26,-1,-1,-1,-1,
//					 -1,-1,-1,-1,-1,-1,-1,-1,-1};
//const double test_problem::stencil[9]={0,0,0,0,64,0,0,0,0};

double test_problem::fix[1]={-32/12.};
