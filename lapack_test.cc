#include <cstdio>
#include <cstdlib>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
	int dgetrs_(char *trans_,int *n,int *nrhs,double *a,int *lda,
		    int *ipiv,double *b,int *ldb,int *info);
	int dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
}

// Solves the matrix system Ax=b, returning the answer in the x array
void solve_matrix(int n,double *A,double *x) {

	// Create the temporary memory that LAPACK needs
	int info,nrhs=1,*ipiv=new int[n];
	char trans='N';

	// Perform the LU decomposition
	dgetrf_(&n,&n,A,&n,ipiv,&info);
	if(info!=0) {
		fputs("LAPACK LU routine failed\n",stderr);
		exit(1);
	}

	// Use the LU decomposition to solve the system
	dgetrs_(&trans,&n,&nrhs,A,&n,ipiv,x,&n,&info);
	if(info!=0) {
		fputs("LAPACK solve routine failed\n",stderr);
		exit(1);
	}

	// Remove temporary memory
	delete [] ipiv;
}

int main() {
	double A[4]={1,1,2,1};
	double b[2]={3,7};

	// Print the matrix and the right hand side
	printf("A=[%6g %6g ]\n  [%6g %6g ]\n\n",*A,A[2],A[1],A[3]);
	printf("b=[%6g ]\n  [%6g ]\n\n",*b,b[1]);

	// Call routine to solve the matrix
	solve_matrix(2,A,b);

	// Print the solution
	printf("x=[%6g ]\n  [%6g ]\n",*b,b[1]);
}
