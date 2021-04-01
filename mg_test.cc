#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "mpi.h"

#include "problem_simple.hh"
#include "multigrid.hh"

// The grid dimensions
const int m=65;
const int n=65;
const int o=65;

int main (int argc, char* argv[]) {
	MPI_Init(&argc,&argv);
	int i;
	comm_buffer com;

	// Set up the processor geometry and problem class
	geometry gm(m,n,o,false,false,false);
	problem_simple p(gm);

	// Set up the region class
	multigrid mg(p,gm,com);

	// Set up the matrices and fields
	mg.setup_matrices(p);
	mg.setup_fields(p);

	// Perform test operations
	mg.setup_test(1,true,3);
	mg.interpolation(3);

	// Make the output directory if it doesn't already exist
	const char odir[]="mgt.odr";
	if(gm.rank==0) mkdir(odir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Output the cross-sections in z
	char buf[64];
	for(int j=2;j<5;j++) {
		if(j<mg.nlevels) {
			for(i=0;i<mg.reg[j]->o;i++) {
				sprintf(buf,"%s/x%d.%d",odir,j,i);
				mg.output_x(buf,i,j);
			}
		}
	}

	// Free the MPI communicators that are part of the geometry and
	// multigrid classes prior to calling MPI_Finalize
	mg.free();gm.free();
	MPI_Finalize();
}
