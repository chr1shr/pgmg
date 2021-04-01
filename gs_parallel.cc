#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "mpi.h"

#include "problem_simple.hh"
#include "buffer.hh"
#include "region.hh"

const int m=120;
const int n=120;
const int o=120;

int main (int argc, char* argv[]) {
	MPI_Init(&argc,&argv);
	int i;
	double err;
	comm_buffer com;

	geometry gm(m,n,o,false,false,false);
	problem_simple p(gm);

	// Set up the region class
	region r(p,gm,com);

	// Set up the matrices and fields
	r.setup_matrices(p);
	r.setup_fields(p);

	// Carry out 200 Gauss--Seidel iterations
	err=r.l2_error();

	double t0=MPI_Wtime();
	if(gm.rank==0) printf("0 %g\n",err);
	for(int i=1;i<=100;i++) {
		r.gauss_seidel();
		err=r.l2_error();
		if(gm.rank==0) printf("%d %g\n",i,err);
	}
	if(gm.rank==0) printf("Time: %g s\n",MPI_Wtime()-t0);

	// Make the output directory if it doesn't already exist
	const char odir[]="gst.odr";
	if(gm.rank==0) mkdir(odir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Output the cross-sections in z
	char buf[64];
	for(i=0;i<o;i++) {
		sprintf(buf,"%s/x.%d",odir,i);
		r.output_x(buf,i);
	}

	gm.free();
	MPI_Finalize();
}
