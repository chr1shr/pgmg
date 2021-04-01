#include "mpi.h"
#include <cstdio>
#include "geometry.hh"
#include "common.hh"
#include "multigrid.hh"
#include "test_problem.hh"

int main(int argc,char* argv[]) {

	MPI_Init(&argc,&argv);

	// system size
	int m = 256;

	// sample problem
	comm_buffer cm;
	geometry gm(m,m,m,false,false,false);
	test_problem pr(gm);
	multigrid mg(pr,gm,cm);

	// set up problem
	mg.setup_matrices(pr);

	// central matrix element of problem_simple is 32/12 ~ 3
	mg.guide.set_tol(3);

	double resid;

	double t0 = MPI_Wtime();

	mg.setup_fields(pr);
	int result = mg.solve_v_cycles(resid);

	double t1 = MPI_Wtime();

	if (result != 0) p_fatal_error("no convergence",-1);
	if (gm.rank==0) printf("first solve: %g\n",t1-t0);

	t0 = MPI_Wtime();
	mg.setup_fields(pr);
	result = mg.solve_v_cycles(resid);
	t1 = MPI_Wtime();

	if (result != 0) p_fatal_error("no convergence",-1);
	if (gm.rank==0) printf("second solve: %g\n",t1-t0);

	t0 = MPI_Wtime();
	mg.setup_fields(pr);
	result = mg.solve_v_cycles(resid);
	t1 = MPI_Wtime();

	if (result != 0) p_fatal_error("no convergence",-1);
	if (gm.rank==0) printf("third solve: %g\n",t1-t0);

	mg.free();
	MPI_Finalize();
}
