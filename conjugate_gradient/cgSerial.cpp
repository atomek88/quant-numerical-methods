// Created by Tomasz Michalik on 5/17/18.
// congugate gradient serial implementation
#include "cg_header.h"

int main(int argc, char* argv[]) {
    // double ** A; double* xvec; double* pvec; double* rvec; double* zvec;
    int i,j,k, nprocs, rank, n, N, x, mypart;

  // #ifdef MPI
    MPI_Status status;
    MPI_Comm comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    int  p = sizeof(double);
  //  #endif
    if (argc != 3){ cli_error();}
    else { n = atoi(argv[1]);}
    // static_assert(n > 0 && isfinite(n), "No good, bad input");
    printf("Solving Poisson Equation on %d x %d domain \n", n,n);
    if (n % nprocs != 0) { printf("Not valid dimensions for this decomp scheme\n"); exit(0);}
    if (!strcmp(argv[2], "serial")){
       // run_dense(n);
        // run sparse
        run_sparse(n);
        } else if (!strcmp(argv[2], "parallel")) {
            // run parallel CG solve
            //printf("My rank %d, mypart:%d starting parallel execution with %d ps\n", rank, mypart, nprocs);
            run_parallel_sparse(n, rank, nprocs, status);
    } else cli_error();
#ifdef MPI
    MPI_Finalize();
#endif
    return 0;
}

void run_dense(long n){
    printf("Solving Dense version of CG\n");
    int N = n*n;
    // alloc matrices
    double** A = matrix(N);
    double* bvec = (double*) malloc(N * sizeof(double));
    double* xvec = (double*) malloc(N * sizeof(double));
    printf("Dense Memory = %.2lf MB\n", (N*N+5*N) * sizeof(double) / 1024.0 / 1024.0);
    printf("Dense Memory = %f MB\n", (((N*N)+(5*N)) * sizeof(double) / 1024.0 ) / 1024.0);
    // compute elements of A possion operator
    fill_A(A, N);
    fill_b(bvec, N);
    // Run dense CG solver
    double start = get_time();
    cg_dense(A, xvec, bvec, N);
    double stop = get_time();
    printf("Dense runtime = %.2lf seconds\n", stop-start);
    // save solution vector to file
    print_matrix(A, N);
    save_vector(xvec, N, "dense.out");
    // free A matrix
    matrix_free(A);
    // free vectors
    free(xvec);
    free(bvec);
}
// sparse
void run_sparse(long n) {
    printf("Running sparse solver CG\n");
    int N = n*n;
    // vector flush
    double* xvec = (double*) calloc(N, sizeof(double));
    double* bvec = (double*) calloc(N, sizeof(double));
    printf("Sparse Memory  = %.2lf MB\n", ((5*N)*sizeof(double) /1024.0) /1024.0);
    // compute elements on boundary condition vector b
    fill_b(bvec, N);
   // print_vector(bvec, N); // test print to compare parallel
    //run sparse solver
    double start = get_time();
    cg_sparse_poisson(xvec, bvec, N);
    double stop = get_time();
    printf("Sparse runtime = %.2lf seconds\n", stop-start);
  //  save
    save_vector(bvec, N, "BEGsparse.out");
    save_vector(xvec, N, "sparse.out");
    free(xvec);
    free(bvec);
}
// run parallel solve sparse
void run_parallel_sparse( long n, int rank, int nprocs, MPI_Status status ){ // pass n is already divided by
    // calculate decomposition
    int N = n*n;
    MPI_File file;
    int mypart = N/ nprocs;
    printf("Sparse Parallel version N:%d,  part:%d...\n", N, mypart);
    double* xvec = (double*) calloc(mypart, sizeof(double)); // these need to be smaller
    double* bvec = (double*) calloc(mypart, sizeof(double));
    double start = get_time();
    // fill and compute elements on boundary condition vector b
    parallel_fill_b( bvec, N, rank, nprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    parallel_cg_sparse_poisson(xvec, bvec, N, rank, nprocs, mypart, status);
    double stop = get_time();
    printf("Sparse runtime = %.2lf seconds\n", stop-start);
    // save init conditions
    //save_vector_p(bvec, N, rank, nprocs, "BEGpar_sparse.out");
    //MPI_Barrier(MPI_COMM_WORLD);
    //save_vector_mpi( bvec,  N,  rank, nprocs, file,  status);
    // save end conditions
    save_vector_p(xvec, N, rank, nprocs, "par_sparse.out");
    MPI_Barrier(MPI_COMM_WORLD);
    free(xvec);
    free(bvec);
    MPI_Finalize();
    return;
}