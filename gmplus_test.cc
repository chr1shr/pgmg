#include "mpi.h"
#include "geometry.hh"
#include "problem_simple.hh"
#include "buffer.hh"
#include "multigrid.hh"

const int m = 100;
const int n = m;
const int o = m;
const bool x_prd = false;
const bool y_prd = false;
const bool z_prd = false;

int main(int argc, char* argv[]) {

	// initialize
	MPI_Init(&argc,&argv);
	geometry gm100(m,n,o,x_prd,y_prd,z_prd);
	geometry gm101(m+1,n+1,o+1,x_prd,y_prd,z_prd);
	geometry gm100p(gm100);

	problem_simple p100(gm100);
	problem_simple p101(gm101);
	problem_simple p100p(gm100p);

	comm_buffer cb;

	/*
	for (int i = 0; i < gm_plus.procs; i++) {

		// print node ownership
		if (gm_plus.rank == i) printf("\nproc %02d (orig %02d) owns...\n"
			"    x: %2d->%2d (%2d nodes), originally (%2d->%2d) (%2d nodes)\n"
			"    y: %2d->%2d (%2d nodes), originally (%2d->%2d) (%2d nodes)\n"
			"    z: %2d->%2d (%2d nodes), originally (%2d->%2d) (%2d nodes)\n",
			gm_plus.rank,gm_root.rank,
			gm_plus.ai,gm_plus.bi,gm_plus.sm,gm_root.ai,gm_root.bi,gm_root.sm,
			gm_plus.aj,gm_plus.bj,gm_plus.sn,gm_root.aj,gm_root.bj,gm_root.sn,
			gm_plus.ak,gm_plus.bk,gm_plus.so,gm_root.ak,gm_root.bk,gm_root.so);

		// sync up so it's readable
		MPI_Barrier(world);
	}
	*/

	if (gm100.rank == 0)
		printf("building multigrid using native (%d,%d,%d) geometry and problem...\n",
			m,n,o);
	multigrid mg100(p100,gm100,cb);
	MPI_Barrier(world);

	if (gm100.rank == 0)
		printf("building multigrid using native (%d,%d,%d) geometry and problem...\n",
			m+1,n+1,o+1);
	multigrid mg101(p101,gm101,cb);
	MPI_Barrier(world);

	if (gm100.rank == 0)
		printf("bulding multigrid using \"+1\" (%d,%d,%d) geometry and problem\n",m+1,n+1,o+1);
	multigrid mg100p_g(p100p,gm100p,cb);
	MPI_Barrier(world);

	if (gm100.rank == 0)
		printf("building multigrid using \"+1\" (%d,%d,%d) geometry and native (%d,%d,%d) problem\n",m+1,n+1,o+1,m,n,o);
	multigrid mg100p_gp(p100,gm100p,cb);
	MPI_Barrier(world);

	gm100.free();
	gm101.free();
	gm100p.free();
	MPI_Finalize();
}
