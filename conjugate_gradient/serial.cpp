//
// Created by Tomasz Michalik on 5/17/18.
//

#include "cg_header.h"

void cg_dense(double ** A, double * x, double * b, long N)
{
    // r = -A*x + b
    double * r = (double *) malloc( N*sizeof(double));
    matvec(r, A, x, N);
    axpy(-1.0, r, 1.0, b, N);

    //p = r;
    double * p = (double *) malloc( N*sizeof(double));
    memcpy(p, r, N*sizeof(double));

    //rsold = r' * r;
    double rsold = dotp(r, r, N);

    // z
    double * z = (double *) malloc( N*sizeof(double));

    long iter = 0;
    for( iter = 0; iter < N; iter++ )
    {
        //z = A * p;
        matvec(z, A, p, N);

        //alpha = rsold / (p' * z);
        double alpha = rsold / dotp(p, z, N);

        //x = x + alpha * p;
        axpy(1.0, x, alpha, p, N);

        //r = r - alpha * z;
        axpy(1.0, r, -alpha, z, N);

        double rsnew = dotp(r,r,N);

        if( sqrt(rsnew) < 1.0e-10 )
            break;

        //p = (rsnew / rsold) * p + r;
        axpy(rsnew/rsold, p, 1.0, r, N);

        rsold = rsnew;
    }
    printf("CG converged in %ld iterations.\n", iter);

    free(r);
    free(p);
    free(z);
}

// Serial Conjugate Gradient Solver Function for Ax = b
// where A is the 2D poisson matrix operator that
// is assembled on the fly (not stored)
// x = 1D solution vector
// b = 1D vector
// N = dimension
void cg_sparse_poisson(double * x, double * b, long N)
{
    // r = -A*x + b
    double * r = (double *) malloc( N*sizeof(double));
    matvec_OTF(r, x, N);
    axpy(-1.0, r, 1.0, b, N);

    //p = r;
    double * p = (double *) malloc( N*sizeof(double));
    memcpy(p, r, N*sizeof(double));

    //rsold = r' * r;
    double rsold = dotp(r, r, N);

    printf("rsold:%.9le\n",rsold);
    // z
    double * z = (double *) malloc( N*sizeof(double));

    long iter = 0;
    for( iter = 0; iter < N; iter++ )
    {
        //z = A * p;
        matvec_OTF(z, p, N);

        //alpha = rsold / (p' * z);
        double t = dotp(p, z, N);
        double alpha = rsold / t;
        //x = x + alpha * p;
        axpy(1.0, x, alpha, p, N);

        //r = r - alpha * z;
        axpy(1.0, r, -alpha, z, N);

        double rsnew = dotp(r,r,N);

        printf("count[%d] sum:%.9le alpha:%.9le rsnew:%.9le\n",  iter,t,alpha, rsnew);
        if( sqrt(rsnew) < 1.0e-10 )
            break;

        //p = (rsnew / rsold) * p + r;
        axpy(rsnew/rsold, p, 1.0, r, N);

        rsold = rsnew;
        printf("rsold is now:%.9le\n",rsold);

    }
    printf("CG converged in %ld iterations.\n", iter);

    free(r);
    free(p);
    free(z);
}

// General Matrix Vector Product for v = M * w
// v, w = 1D vectors
// M = 2D matrix
// N = dimension
void matvec( double * v, double ** M, double * w, long N )
{
    // Set solution vector to 0
    memset( v, 0, N*sizeof(double));

    for (long i = 0; i < N; i++)
        for (long j = 0; j < N; j++)
            v[i] += (M[i][j] * w[j]); // makes sense of multiplying each matrix indexed
}

// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void matvec_OTF( double * v, double * w, long N )
{   // exchange information among neighbors
    // Determine physical dimension
    // send each neighbor a n length vector at start
    long n = sqrt(N);

    // Set solution vector to 0
    memset( v, 0, N*sizeof(double));

    for( long i = 0; i < N; i++ )
    {   // send data again to neighbors/ bcast?
        // Physical 2D x-coordinate
        long x_i = i%n;
        printf("i:%d, x_i:%d my first el: %.9le \n", i, x_i, w[i]);
        // Far left diagonal
        double far_left_diag = 0.0;
        if( i >= n ) {
            far_left_diag = -1.0 * w[i - n];
          //  printf("after far left  fld%.9le  w[i-n]:%.9le \n", far_left_diag, w[i - n]);
        }
        // left diagonal
        double left_diag = 0.0;
        if( x_i != 0 ) {
            left_diag = -1.0 * w[i - 1];
      //      printf("after left  ld:%.9le  w[i-1]:%.9le \n", far_left_diag, w[i-1]);
        }
        // Main diagonal
        double main_diag = 4.0 * w[i];

        // Right diagonal
        double right_diag = 0.0;
        if( x_i != n-1 ) {
            right_diag = -1.0 * w[i + 1];
         //   printf("after right  rd%.9le  w[i+1]:%.9le \n", right_diag, w[i+1]);
        }
        // Far right diagonal
        double far_right_diag = 0.0;
        if( i < N - n ) {
            far_right_diag = -1.0 * w[i + n];
      //      printf("after far right frd:%.9le  w[i+n]:%.9le \n", far_right_diag, w[i+n]);
        }

        v[i] = far_left_diag + left_diag + main_diag + right_diag + far_right_diag;
        printf("v[%d]:%.9le \n",i, v[i]);
        printf("farleft:%.9le left:%.9le mid:%.9le right:%.9le farright:%.9le \n",far_left_diag,left_diag,main_diag,right_diag,far_right_diag);
    }
}

// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double dotp( double * a, double * b, long N) {
    double c = 0.0;
    for( long i = 0; i < N; i++ ){
        c += a[i]*b[i];
     //   printf("dotp i:%d a:%.9le b:%.9le sum:%.9le\n",i,c);
        }

    return c;
}

// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void axpy( double alpha, double * w, double beta, double * v, long N) {
    for( long i = 0; i < N; i++ ) {
        w[i] = alpha * w[i] + beta * v[i];
       // if(i %n ==0){
        printf("axpy,i:%d beta:%.9le w[i]:%.9le v[i]:%.9le \n",i, beta, i, w[i], v[i]);
           //}

    }
}


// Fills an Explicit Fully Stored 2D Poisson Operator Matrix 'A'
// A = 2D Matrix
// N = dimension
void fill_A(double ** A, long N)
{
    long n = sqrt(N);

    for( long i = 0; i < N; i++ )
    {
        for( long j = 0; j < N; j++ )
        {
            // Compute x-coordinate of index
            long x = j % n;

            // main diagonal
            if(i == j)
                A[i][j] = 4.0;
                // left diagonal
            else if (i == j+1 )
            {
                A[i][j] = -1.0;
                if(x == n-1)
                    A[i][j] = 0.0;
            }
                // right diagonal
            else if (i == j-1 )
            {
                A[i][j] = -1.0;
                if(x == 0)
                    A[i][j] = 0.0;
            }
                // far left diagonal
            else if(j+n == i)
                A[i][j] = -1.0;
                // far right diagonal
            else if(j-n == i)
                A[i][j] = -1.0;
                // Otherwise, A is 0
            else
                A[i][j] = 0.0;
        }
    }
}
// Sets a circular source term for the right hand side 'b' vector
// Takes as input 2D spatial domain indices
// i, j = indices
// n = physical dimension
double find_b(long i, long j, long n)
{  // printf("In serial find/get b, i:%d, j:%d n:%d \n", i,j,n);
    double delta = 1.0 / (double) n;

    double x = -.5 + delta + delta * j;
    double y = -.5 + delta + delta * i;
    // Check if within a circle
    double radius = 0.1;
  //  printf("x: %2.1lf y:%2.1lf delta:%2.1lf \n", x,y,delta);
    if( x*x + y*y < radius*radius )
        return delta * delta / 1.075271758e-02;
    else
        return 0.0;
}
// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void fill_b(double * b, long N)
{  // printf("size of N:%d \n", N);
    long n = sqrt(N);
    for( long i = 0; i < n; i++ )
        for( long j = 0; j < n; j++ )
        {
            long idx = i*n + j;
            b[idx] = find_b(i,j,n);
        //    printf("i:%d j:%d idx:%d b:%.9le \n",i,j,idx, b[idx]);
        }
   // printf("\n");

}