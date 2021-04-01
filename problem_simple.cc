/** \file problem_simple.cc
 * \brief Function implementations for the problem_simple class. */

#include "common.hh"
#include "problem_simple.hh"

/** Fills the entries in the table with stencil values for a particular gridpoint.
 * \param[in] (i,j,k) the index of the gridpoint to consider.
 * \param[in,out] en a reference to the pointer to the stencil table. */
void problem_simple::fill_entries(int i,int j,int k,double *&en) {
	p_fatal_error("No need to fill entries with this problem type",1);
}

/** An array containing the basic finite-difference stencil for the Laplacian
 * operator. */
//double problem_simple::stencil[27]={-1,-1,-1,-1,0,-1,-1,-1,-1,-1,0,-1,0,6,0,-1,0,-1,-1,-1,-1,-1,0,-1,-1,-1,-1};
// Nodal Poisson problem
// double problem_simple::stencil[27]={-1,-2,-1,-2,0,-2,-1,-2,-1,-2,0,-2,0,32,0,-2,0,-2,-1,-2,-1,-2,0,-2,-1,-2,-1};
//double problem_simple::stencil[27]={1/12.,2/12.,1/12.,2/12.,0,2/12.,1/12.,2/12.,1/12.,2/12.,0,2/12.,0,-32/12.,0,2/12.,0,2/12.,1/12.,2/12.,1/12.,2/12.,0,2/12.,1/12.,2/12.,1/12.};
// Cell center Poisson problem
double problem_simple::stencil[27]={0,0,0,0,1,0,0,0,0,0,1,0,1,-6,1,0,1,0,0,0,0,0,1,0,0,0,0};
//double problem_simple::stencil[27]={0,0,0,0,0,0,0,0,0,
//					0,0,0,1,-2,1,0,0,0,
//					0,0,0,0,0,0,0,0,0};
//const double problem_simple::stencil[27]={1,-2,1,-2,4,-2,1,-2,1,
//					-2,4,-2,4,-8,4,-2,4,-2,
//					1,-2,1,-2,4,-2,1,-2,1};
//const double problem_simple::stencil[27]={-1,-1,-1,-1,-1,-1,-1,-1,-1,
//					 -1,-1,-1,-1,26,-1,-1,-1,-1,
//					 -1,-1,-1,-1,-1,-1,-1,-1,-1};
//const double problem_simple::stencil[9]={0,0,0,0,64,0,0,0,0};

double problem_simple::fix[1]={-0.375};
