//
// Created by Tomasz Michalik on 5/17/18.
//

#include "cg_header.h"
// Parallel MPI version conjugate gradient solver Ax = b
// A = 2d posson matrix operator // mypart = size for each proc/rank
// x = 1d solution , b = 1dvector
void parallel_cg_sparse_poisson(double* x, double* b, long N, int rank, int nprocs, int mypart, MPI_Status status){
   // printf("My rank %d, in function par sparse poisson for %d with %d procs partition:%d\n", rank, N, nprocs,mypart);
    double* nbuf = (double*) calloc(sqrt(N), sizeof(double));
    double* sbuf = (double*) calloc(sqrt(N), sizeof(double)); // load 0..n and last n elements
    // r = -A*x + b
    double* r = (double*) malloc (mypart*sizeof(double));
    parallel_matvec_OTF(r, x, N, rank, nprocs, status, nbuf, sbuf); //- complicated ghost exchange
    // parallel axpy
   // MPI_Barrier(MPI_COMM_WORLD);
    parallel_axpy(-1.0, r, 1.0, b, N, rank, nprocs, -1);
    // p = r
    double* p = (double*) malloc (mypart*sizeof(double));
    memcpy(p, r, mypart*sizeof(double));
    // r' * r dotp
    double rsold = parallel_dotp(r, r, N, rank, nprocs, -1);
  //  MPI_Barrier(MPI_COMM_WORLD);
    // z
    double* z = (double*) malloc (mypart*sizeof(double));
    long count = 0;
    for( count = 0; count< N; count++) { /// change count to N, 4 for testing
        // z = A * p
        parallel_matvec_OTF(z, p, N, rank, nprocs, status, nbuf, sbuf);
        // alpha = rsold / p' * z
       // MPI_Barrier(MPI_COMM_WORLD);
        double oldalpha = parallel_dotp(p, z, N, rank, nprocs, count);
        double alpha = rsold / oldalpha;
      //  MPI_Barrier(MPI_COMM_WORLD);
        // x = x+ alpha * p
        parallel_axpy(1.0,x, alpha, p, N, rank, nprocs, count);
        // r = r- alpha * z;
        parallel_axpy(1.0, r, -alpha, z, N, rank, nprocs, count);
       // MPI_Barrier(MPI_COMM_WORLD);
        double rsnew = parallel_dotp(r, r, N, rank, nprocs, count);
       // MPI_Barrier(MPI_COMM_WORLD);
        if (sqrt(rsnew) < 1.0e-10 ) {
            printf("break cond: %.9le\n", sqrt(rsnew));
            break; // break conidtion
        }
        // p = rsnew / rsold * p + r
        parallel_axpy(rsnew/rsold, p, 1.0, r, N, rank, nprocs, count);
        rsold = rsnew;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    printf("CG converged in %ld iterations\n", count);
    free(r); free(p); free(z);
}
// Parallel Distributed MPI Version of
// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void parallel_matvec_OTF( double * v, double * w, long N, int rank, int nprocs, MPI_Status status, double* nbuf, double* sbuf )
{   // physical dimension, but send n size vector to my neighbor
    long n = sqrt(N);
    int mypart = N/nprocs; int i,j; // my part of matvec
    int offset = mypart - n;
    // printf("in parallel matvec N:%d rank%d n:%d mypart:%d \n", N, rank, n, mypart);
    // buffers for north and south, n sized

    // solution vector set to 0
    memset (v, 0, mypart*sizeof(double));
    int north = rank - 1; if (north < 0) north = -1; // MPI_PROC_NULL
    int south = rank + 1; if (south >=nprocs) south = -1;
    // printf("My rank:%d N:%d n:%d my north: %d my south:%d length:%d\n", rank, N, n, north, south, offset);
// i only need to exchange once at beginning of routine to send each others neighbors. need to update v?
    if (rank == 0) {
        //  printf("My rank:%d sending south, receiving south offset:%d \n", rank, offset);
        MPI_Sendrecv(&w[offset], n, MPI_DOUBLE, south, 11, sbuf, n, MPI_DOUBLE, south, 4, MPI_COMM_WORLD, &status);
    }
    if (rank >0 && rank < nprocs-1) {
        MPI_Sendrecv(&w[0], n, MPI_DOUBLE, north, 4, nbuf, n, MPI_DOUBLE, north, 11, MPI_COMM_WORLD, &status);
        // printf("My rank:%d sending both, receiving both., val:%f with amp:%f next:%f \n", rank, w, &w[0], w[1]);
        MPI_Sendrecv(&w[offset], n, MPI_DOUBLE, south, 11, sbuf, n, MPI_DOUBLE, south, 4, MPI_COMM_WORLD, &status);
    }
    if (rank == nprocs - 1) {

        MPI_Sendrecv(&w[0], n, MPI_DOUBLE, north, 4, nbuf, n, MPI_DOUBLE, north, 11, MPI_COMM_WORLD, &status);
        // printf("My rank:%d sending north, receiving northval:%.9le with amp:%f next:%f \n", rank, w[1], nbuf[0], nbuf[1]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for ( long i=0; i<mypart; i++) { // use N or mypart here?
        // physical 2d xcoord

        long x_i = i % n;
        // Far left diagonal
        long my_i = i + (rank*mypart); long tp = my_i  % n;
        //   printf("me%d i:%d, x_i:%d my_i:%d my first el: %.9le tp:%d \n",rank, i, x_i, my_i, w[i], tp);
        // Far left diagonal
        double far_left_diag = 0.0;
        if( my_i >= n ) {
            if (i < n){
                //   printf("took nbuf\n");
                far_left_diag = -1.0 * nbuf[tp];
            } else {
                //      printf("took w-i\n");
                far_left_diag = -1.0 * w[i-n];
            }

            //  printf("after far left  fld:%.9le rank%d nbuf[%d]:%.9le i-n:%.9le \n", far_left_diag, rank,tp, nbuf[x_i], w[i-n]);
        }
        // left diagonal
        double left_diag = 0.0;
        // alternate
        // if (my_i != 0) { //( x_i  == 0 && my_i != 0){
        if (tp != 0 ){
            left_diag = -1.0 * w[i - 1];
        }
        // Main diagonal
        double main_diag = 4.0 * w[i];

        // Right diagonal, keep contiguous except for ends
        double right_diag = 0.0;
        //  if (my_i  != N ){ // wraparound for contiguity, can % mypart?
        if (tp != n-1) {
            right_diag = -1.0 * w[i + 1];
            //          printf("after r1 rightd  rd:%.9le rank%d w[i+1]:%.9le \n", right_diag, rank, w[i + 1]);
        }
        // Far right diagonal
        double far_right_diag = 0.0;
        if( my_i < N - n ) {
            if (i < mypart - n) { // use only if not on southern border
                far_right_diag = -1.0 * w[i + n];
                //    printf("after1 far rightd myx_i:%d frd:%.9le rank%d w[i_n[%d]:%.9le \n", x_i, far_right_diag, rank, tp,w[i+n]);

            } else {
                far_right_diag = -1.0 * sbuf[x_i];
                //    printf("after2 far rightd myx_i:%d frd:%.9le rank%d sbuf[%d]:%.9le \n", x_i, far_right_diag, rank, tp,sbuf[x_i]);
            }

        }
       // MPI_Barrier(MPI_COMM_WORLD);
        v[i] = far_left_diag + left_diag + main_diag + right_diag + far_right_diag;
       // MPI_Barrier(MPI_COMM_WORLD);


    }

}

// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension, dimension for MYPART???
double parallel_dotp( double * a, double * b, long N, int mype, int nprocs, int iter) // using rank and nprocs compute where your dotp is
{
    int mypart = N/nprocs;
    long i; double localsum = 0.0; double finsum = 0.0;

    for (i=0; i<mypart; i++){ // use N here or mypart?
        localsum += a[i] * b[i];

    }
  //  MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&localsum, &finsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   // MPI_Barrier(MPI_COMM_WORLD);

    return finsum;
}

// Parallel Distributed MPI Version of
// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void parallel_axpy( double alpha, double * w, double beta, double * v, long N, int mype, int nprocs, int iter) {
    int  n = sqrt(N);
    int mypart = N / nprocs;

    for (long i=0; i< mypart; i++){
        w[i] = (alpha * w[i]) + (beta * v[i]);

    }
}
// Source term find_b method
double par_find_b(long i,long j, long n){
    //printf("In parallel find/get b, i:%d, j:%d n:%d \n", i,j,n);
    double delta = 1.0 / (double) n;
    double x = -.5 + delta + delta * j;
    double y = -.5 + delta + delta * i;
    double radius = 0.1;
    //  printf("x: %2.1lf y:%2.1lf delta:%2.1lf \n", x,y,delta);
    if (x*x + y*y < radius*radius) {
        return delta*delta /1.075271758e-02;
    }
    else {return 0.0;}
}
// Parallel Distributed MPI Version of
// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector      // each proc has N/proc vector length
// N = dimension (length of b)
void parallel_fill_b(double * b, long N, int mype, int nprocs)
{
    long n = sqrt(N);
    int mypart = N / nprocs; // each part of each proc
    int size = n/nprocs;
    //printf("** par fill b, me:%d, N:%d,part:%d \n", mype, N,  mypart );
    for (long i = 0; i<size; i++){
        for (long j = 0; j< n; j++) {
            long idx = i *n +j ;
            long my_x = i+ (mype * size);
            b[idx] = par_find_b(my_x, j, n);
            // printf("rank:%d i:%d j:%d myx:%d idx:%d b:%.9le \n", i, j, mype,my_x, idx, b[idx]);
        }


    }

}
