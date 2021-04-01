#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <cmath>
#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <iostream>
#include "mpi.h"
#include "problem_simple.hh"
#include "multigrid.hh"
/** A general class of linear operator y = A * x.
 * Spectra needs rows(), cols() and perform_op()
 */
class LinearOp {

    public:
        problem_simple * p;
        const int m;
        const int n;
        const int o;
        const int mn;
        const int N;

        LinearOp (problem_simple *p_) : p(p_),
        m(p->m), n(p->n), o(p->o), mn(m*n),
        N(p->m * p->n * p->o) {}

        // return index in the vector that corresponds to grid point (i,j,k)
        int index(int i, int j, int k) {
            return (i%m)+(j%n)*m+(k%o)*mn;
        }
        // return index in the 27 vector that corresponds to stencil entry coordinate
        // [-1,1] x [-1,1] x [-1,1]
        int st_index(int si, int sj, int sk) {
            return (si+1) + 3*(sj+1) + 9*(sk+1);
        }

        /** The following 3 functions are needed by Spectra */
        virtual int rows() {return N;}
        virtual int cols() {return N;}
        // y_out = M * x_in
        virtual void perform_op(const double *x_in, double *y_out) {
            for (int k=0;k<o;k++) {
                for(int j=0;j<n;j++) {
                    for(int i=0;i<m;i++){
                        int ind = index(i,j,k);
                        if(p->interior(i,j,k)) {
                            y_out[ind] = 0.;
            // Loop over the stencil
            for (int sk=-1;sk<=1;sk++) {
                for (int sj=-1;sj<=1;sj++) {
                    for (int si=-1;si<=1;si++){
                        y_out[ind] += p->stencil[st_index(si,sj,sk)] * x_in[index(i+si, j+sj, k+sk)];
                    }
                }
            }
                        } else {
                            y_out[ind] = p->fix[0] * x_in[ind];
                        }
                    }
                }
            }
        }
};

/** A class that uses PGMG to define the linear operation of y = A * x.
 * Since the x_in and y_out used by Spectra doesn't have ghost cells/nodes, we have to define
 * a new function in PGMG, as well as region class, that does A*x without the padded data.
 * Also note that this can only work on 1 MPI process. Periodic condition is taken care of
 * by using a non-negative mod function to wrap the indices around.*/
class MGLinearOp : public LinearOp {
    public:
        const int lev;
        multigrid *mg;
        MGLinearOp(problem_simple *p_, multigrid * mg_, const int lev_) : LinearOp(p_), lev(lev_), mg(mg_){}

        virtual int rows () {
            return mg->reg[lev]->m* mg->reg[lev]->n*mg->reg[lev]->o;
        }
        virtual int cols () {
            return mg->reg[lev]->m* mg->reg[lev]->n*mg->reg[lev]->o;
        }
        virtual void perform_op(const double *x_in, double *y_out) {
           mg->unpadded_Ax(x_in, y_out, lev);
        }
};

/* Tell the compiler about LAPACK functions; here we use dense symmetric eigenvalues solver.*/
extern "C" {
    void dsyev_(
            char *JOBZ,
            char *UPLO,
            int 	*N,
            double  *A,
            int 	*LDA,
            double  *W,
            double *WORK,
            int *LWORK,
            int* INFO
   );

    void dgeev_(
            char *JOBVL,
            char *JOBVR,
            int 	*N,
            double  *A,
            int 	*LDA,
            double  *WR,
            double  *WL,
            double  *VL,
            int     *LDVL,
            double  *VR,
            int     *LDVR,
            double  *WORK,
            int     *LWORK,
            int     *INFO
   );
}

/** Computes the eigenvalues, and optionally eigenvectors, of the n by n
 * symmetric matrix A, storing the results in w. The lower triangular
 * part of the input matrix is ignored. */
void eigs_symmetric(MGLinearOp *mg, int n, bool evecs, const char* filename) {
    int info, ldab=n, lwork=3*n;
    char uplo='U';
    char jobz;
    if(evecs) jobz='V';
    else jobz='N';

    // Reorganize the array into the correct format
    double *A =new double[n*n];
    double *n1d_in = new double[n];
    double *n1d_out = new double[n];
    double *ab=new double[ldab*n];
    double *work=new double[lwork];
    for(int i=0;i<n*n;i++) ab[i] = 0.;
    for(int i=0;i<n;i++) {
            n1d_in[i] = 0.;
            n1d_out[i] = 0.;
    }
    // We use perform_op function in the linear operator to populate A
    for (int col=0;col<n;col++) {
        n1d_in[col] = 1.0;
        // set the previous element back to zero
        if(col>0) n1d_in[col-1] = 0.;
        mg->perform_op(n1d_in, n1d_out);
        // For symmetric matrices, row and column can be switched
        // But we keep them ordered, in case we adapt this for general matrices
        for (int row=0;row<n;row++){
        //   printf("col %d row %d n1d_out %g\n", col, row, n1d_out[row]);
            A[col+row*n] = n1d_out[row];
        }
    }

    /*
    //Print out the entire matrix for inspection
    if(n<220) {
        for(int row=0;row<n;row++){
            for(int col=0;col<n;col++){
                printf("%g ", A[row*n + col]);
            }
            printf("\n");
        }
    }
    return ;
    */

    // notice that (j,k) is (row, column) but lapack store the matrix column wise
    for(int k=0;k<n;k++) {
        for(int j=0;j<=k;j++) {
            // For symmetric matrix, just the upper or the lower triangular matrix
            ab[j+k*n] = A[k+j*n];
        }
    }

    // zero out the output array
    for(int i=0;i<n;i++) {
            n1d_out[i] = 0.;
    }

    // Make call to LAPACK routine
    dsyev_(&jobz,&uplo,&n,ab,&ldab,n1d_out,work,&lwork,&info);
    if(info!=0) {
        fputs("LAPACK routine failed\n",stderr);
        exit(1);
    }

    // Do something with the eigenvalues
    // Maybe we should sort them by magnitude
    FILE *fh = safe_fopen(filename, "w");
    fprintf(fh, "# System size %d\n", n);
    double max=-1e100, min=1e100;
    int maxsign=1, minsign=1;
    for (int i=0;i<n;i++) {
        double mag = fabs(n1d_out[i]);
        fprintf(fh, "%12.8g\n", n1d_out[i]);
        if(mag > max) {max= mag; maxsign=int(n1d_out[i]/mag);}
        if(mag < min) {min= mag; minsign=int(n1d_out[i]/mag);}
    }
    fclose(fh);

    printf("Done. See file %s.\n"
           "Sneak peak: max (%g) and min (%g) of the eigenvalues\n", filename, max*maxsign, min*minsign);

    // Remove temporary memory
    delete [] A;
    delete [] ab;
    delete [] work;
    delete [] n1d_in;
    delete [] n1d_out;
}

/** Computes the eigenvalues, and optionally eigenvectors, of the n by n
 * symmetric matrix A, storing the results in w. The lower triangular
 * part of the input matrix is ignored. */
void eigs_gen(MGLinearOp *mg, int n, const char* filename) {
    int info, ldab=n, lwork=3*n;
    char jobvl='N';
    char jobvr='N';
    int ldvl = 1, ldvr = 1;
    double *fake_v_array = new double[1];

    // Reorganize the array into the correct format
    double *A =new double[n*n];
    // The folllowing two array will double as real and imaginary part of the eigenvalues
    double *n1d_in = new double[n];
    double *n1d_out = new double[n];
    double *ab=new double[ldab*n];
    double *work=new double[lwork];
    for(int i=0;i<n*n;i++) ab[i] = 0.;
    for(int i=0;i<n;i++) {
            n1d_in[i] = 0.;
            n1d_out[i] = 0.;
    }
    // We use perform_op function in the linear operator to populate A
    for (int col=0;col<n;col++) {
        n1d_in[col] = 1.0;
        // set the previous element back to zero
        if(col>0) n1d_in[col-1] = 0.;
        mg->perform_op(n1d_in, n1d_out);
        // For symmetric matrices, row and column can be switched
        // But we keep them ordered, in case we adapt this for general matrices
        for (int row=0;row<n;row++){
        //   printf("col %d row %d n1d_out %g\n", col, row, n1d_out[row]);
            A[col+row*n] = n1d_out[row];
        }
    }

    /*
    //Print out the entire matrix for inspection
    if(n<220) {
        for(int row=0;row<n;row++){
            for(int col=0;col<n;col++){
                printf("%g ", A[row*n + col]);
            }
            printf("\n");
        }
    }
    return ;
    */

    // notice that (j,k) is (row, column) but lapack store the matrix column wise
    for(int k=0;k<n;k++) {
        for(int j=0;j<n;j++) {
            // For asymmetric general matrix, we populate the whole matrix
            ab[j+k*n] = A[k+j*n];
        }
    }

    // zero out the output array
    for(int i=0;i<n;i++) {
            n1d_in[i] = 0.;
            n1d_out[i] = 0.;
    }

    // Make call to LAPACK routine
    dgeev_(&jobvl, &jobvr, &n, ab, &ldab, n1d_in, n1d_out, fake_v_array, &ldvl, fake_v_array, &ldvr, work, &lwork, &info);
    if(info!=0) {
        fputs("LAPACK routine failed\n",stderr);
        exit(1);
    }

    // Do something with the eigenvalues
    // Maybe we should sort them by magnitude
    FILE *fh = safe_fopen(filename, "w");
    fprintf(fh, "# System size %d\n", n);
    fprintf(fh, "# real imaginary\n");
    double max=-1e100, min=1e100;
    int maxsign=1, minsign=1;
    for (int i=0;i<n;i++) {
        double mag = sqrt(n1d_in[i]*n1d_in[i] + n1d_out[i]*n1d_out[i]);
        fprintf(fh, "%12.8g %12.8g\n", n1d_in[i], n1d_out[i]);
        if(mag > max) {max= mag; maxsign=i;}
        if(mag < min) {min= mag; minsign=i;}
    }
    fclose(fh);

    printf("Done. See file %s.\n"
           "Sneak peak: max mag=%g, (%g+%gi) and min mag=%g (%g+%gi) of the eigenvalues\n", filename, max, n1d_in[maxsign], n1d_out[maxsign],
            min, n1d_in[minsign], n1d_out[minsign]);

    // Remove temporary memory
    delete [] A;
    delete [] ab;
    delete [] work;
    delete [] n1d_in;
    delete [] n1d_out;
    delete [] fake_v_array;
}


int main (int argc, char* argv[]) {
	MPI_Init(&argc,&argv);
    // Shared variables, system size m x m x m
    int m = 6;
    int nconv = 0;
    int nev = 3;
    if(argc==2) {
        m = atof(argv[1]);
    }
    comm_buffer com;
    char filename[256];

    // Create a PGMG solver with periodic boundary condition
    geometry gm_prd(m,m,m,true,true,true);
	problem_simple p_prd(gm_prd);
    multigrid mg_prd(p_prd, gm_prd, com);
    mg_prd.setup_matrices(p_prd);

    // For each level, we create a linear operator on that level, then use Spectra's eigenvalue solver to find
    // the largest and the smallest eigen values.
    puts("####### PERIODIC PGMG Solver ##########");
    for (int lev=0;lev<mg_prd.nlevels;lev++){
        MGLinearOp mgop_prd(&p_prd, &mg_prd, lev);
        std::cout << "\nLevel " << lev << " : size of the matrix is " << mgop_prd.rows() << std::endl;

        puts("####### SPECTRA ##########");
        // Create a symmetric eigenvalues solver, and get the top nev eigenvalues with the largest magnitude
        Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, LinearOp> eigs_lg_mg_prd(&mgop_prd, nev, 2*nev);
        eigs_lg_mg_prd.init();
        nconv = eigs_lg_mg_prd.compute();

        // Retrieve prd results
        Eigen::VectorXd evalues_lg_mg_prd;
        if(eigs_lg_mg_prd.info() == Spectra::SUCCESSFUL) {
            evalues_lg_mg_prd = eigs_lg_mg_prd.eigenvalues();
            std::cout << "Converged in " << nconv << " iterations. " << std::endl;
            std::cout << "Largest eigenvalues:\n" << evalues_lg_mg_prd << std::endl;
        } else {
            std::cout << "Not successful in " << nconv << " iterations\n";
        }

        // Repeat for the top nev eigenvalues with the smallest magnitude
        Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, LinearOp> eigs_sm_mg_prd(&mgop_prd, nev, 2*nev);
        eigs_sm_mg_prd.init();
        nconv = eigs_sm_mg_prd.compute();

        // Retrieve prd results
        Eigen::VectorXd evalues_sm_mg_prd;
        if(eigs_sm_mg_prd.info() == Spectra::SUCCESSFUL) {
            evalues_sm_mg_prd = eigs_sm_mg_prd.eigenvalues();
            std::cout << "Converged in " << nconv << " iterations. " << std::endl;
            std::cout << "Smallest eigenvalues:\n" << evalues_sm_mg_prd << std::endl;
        } else {
            std::cout << "Not successful in " << nconv << " iterations\n";
        }

        // Using LAPACK to solve for the entire Eigenspectrum
        sprintf(filename, "prd_ev_lev_%d", lev);
        puts("####### LAPACK ##########");
        eigs_symmetric(&mgop_prd, mgop_prd.rows(), false, filename);
    }
    gm_prd.free();

    // Create a PGMG solver with Dirichlet boundary condition, the same stencil as above
    geometry gm_nonprd(m,m,m,false,false,false);
    problem_simple p_nonprd(gm_nonprd);
    multigrid mg_nonprd(p_nonprd, gm_nonprd, com);
    mg_nonprd.setup_matrices(p_nonprd);

    // For each level, we create a linear operator on that level, then use Spectra's eigenvalue solver to find
    // the largest and the smallest eigen values.
    puts("\n####### NON-PERIODIC PGMG Solver ##########");
    for (int lev=0;lev<mg_prd.nlevels;lev++){
        MGLinearOp mgop_nonprd(&p_nonprd, &mg_nonprd, lev);
        std::cout << "\nLevel " << lev << " : size of the matrix is " << mgop_nonprd.rows() << std::endl;

        puts("####### SPECTRA ##########");
        Spectra::GenEigsSolver<double, Spectra::LARGEST_MAGN, LinearOp> eigs_lg_mg_nonprd(&mgop_nonprd, nev, 2*nev);
        eigs_lg_mg_nonprd.init();
        nconv = eigs_lg_mg_nonprd.compute();
        // retrieve nonprd results, it's complex
        Eigen::VectorXcd evalues_lg_mg_nonprd;
        if(eigs_lg_mg_nonprd.info() == Spectra::SUCCESSFUL) {
            evalues_lg_mg_nonprd = eigs_lg_mg_nonprd.eigenvalues();
            std::cout << "Converged in " << nconv << " iterations. " << std::endl;
            std::cout << "Largest eigenvalues:\n" << evalues_lg_mg_nonprd << std::endl;
        } else {
            std::cout << "Not successful in " << nconv << " iterations\n";
        }

        Spectra::GenEigsSolver<double, Spectra::SMALLEST_MAGN, LinearOp> eigs_sm_mg_nonprd(&mgop_nonprd, nev, 2*nev);
        eigs_sm_mg_nonprd.init();
        nconv = eigs_sm_mg_nonprd.compute();
        // retrieve nonprd results
        Eigen::VectorXcd evalues_sm_mg_nonprd;
        if(eigs_sm_mg_nonprd.info() == Spectra::SUCCESSFUL) {
            evalues_sm_mg_nonprd = eigs_sm_mg_nonprd.eigenvalues();
            std::cout << "Converged in " << nconv << " iterations. " << std::endl;
            std::cout << "Smallest eigenvalues:\n" << evalues_sm_mg_nonprd << std::endl;
        } else {
            std::cout << "Not successful in " << nconv << " iterations\n";
        }

        // Using LAPACK to solve for the entire Eigenspectrum
        sprintf(filename, "nonprd_ev_lev_%d", lev);
        puts("####### LAPACK ##########");
        eigs_gen(&mgop_nonprd, mgop_nonprd.rows(), filename);
    }
	gm_nonprd.free();

	MPI_Finalize();
}
