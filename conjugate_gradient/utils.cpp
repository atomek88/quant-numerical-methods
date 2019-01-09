//
// Created by Tomasz Michalik on 5/17/18.
#include "cg_header.h"

void print_matrix(double ** A, long N )
{
    for( long i = 0; i < N; i++ )
    {
        for( long j = 0; j < N; j++ )
        {
            printf("%4.1lf ", A[i][j]);
        }
        printf("\n");
    }
}
void print_vector_p(double * x, long N, int rank, int nprocs){
    printf("print vector parallel loop for rank:%d\n",rank);
    int n = sqrt(N);
    int mypart = n / nprocs;
    long xstart = rank * mypart;
    long xend = (rank+1) * mypart;
   //int matrix_dim = sqrt(total);
   //int y_dim = matrix_dim / nprocs;
    long idx = 0; // not use
    printf("My rank %d, starting idx:%d xstart: %d xend:%d..printing my vector\n", rank, idx, xstart, xend);
    for( long i = 0; i < mypart; i++ )
    {
        for( long j = 0; j < n; j++ )
        {
            printf("%.9le ", x[idx]);
            idx++;
        }
        printf("\n");
    }
}
void print_vector(double * x, long N )
{   printf("Serial print vector\n");
    long n = sqrt(N);
    long idx = 0;
    for( long i = 0; i < n; i++ )
    {
        for( long j = 0; j < n; j++ )
        {
            printf("%.9le ", x[idx]);
            idx++;
        }
        printf("\n");
    }
}
// mpi save vector
void save_vector_mpi(double * x, long N, int rank, int nprocs, MPI_File file, MPI_Status status){
    MPI_Datatype num_as_string; int n = sqrt(N);
    int mypart = N/nprocs; int chunk = n/nprocs;
    int xstart = mypart * rank;
    printf("start row %d, chunk %d, length %d printing\n", rank, chunk, mypart);
    MPI_Datatype localarray;
    MPI_Type_contiguous(8, MPI_DOUBLE, &num_as_string);
    MPI_Type_commit(&num_as_string);
    int globalsizes[2] = {n, n};
    int localsizes [2] = {chunk, n};
    int starts[2]      = {xstart, 0};
    MPI_Type_create_subarray(2, globalsizes, localsizes, starts,
                             MPI_ORDER_C, num_as_string,
                             &localarray);
    MPI_Type_commit(&localarray);

    MPI_File_open(MPI_COMM_WORLD, "mpiprint.txt",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);
    //MPI_Offset displace = rank * length * sizeof(int); // for static decomp can add row
    MPI_File_set_view(file, 0, MPI_DOUBLE, localarray,
                      "native", MPI_INFO_NULL);
/* Change here to print 'Tmain' to print Binary. All above related to char conversion can be commented out */
    MPI_File_write_all(file, x, n*chunk, num_as_string, &status); // diff between write all and write?
    MPI_File_close(&file);
    MPI_Type_free(&localarray);
    MPI_Type_free(&num_as_string);
}
// parallel save vector, mostly use for testing if send all back to master at end
void save_vector_p(double * x, long N, int rank, int nprocs, char * fname ){
    printf("Save vector parallel loop, rank:%d \n",rank);
    FILE * fp = fopen(fname, "a");
    int n = sqrt(N);
    int width = n / nprocs; // needs to be evenly divisible
   // long xstart = rank * width;
    long idx = 0;
    for( long i = 0; i < width; i++ )
    {
        for( long j = 0; j < n; j++ )
        {
            fprintf(fp, "%.9le", x[idx]);
            idx++;
            if( j != n-1 )
                fprintf(fp, " ");
        }
        if( i != n - 1 )
            fprintf(fp, "\n");
    }
    fclose(fp);
}
// each have part of the solution vector, but may send it all back to master to use this
void save_vector(double * x, long N, char * fname )
{   printf("Save vector loop\n");
    FILE * fp = fopen(fname, "w");
    long n = sqrt(N);
    long idx = 0;
    for( long i = 0; i < n; i++ )
    {
        for( long j = 0; j < n; j++ )
        {
            fprintf(fp, "%.9le", x[idx]);
            idx++;
            if( j != n-1 )
                fprintf(fp, " ");
        }
        if( i != n - 1 )
            fprintf(fp, "\n");
    }
    fclose(fp);
}

/* Allocates 2-D Contiguous Matrix */
double ** matrix( long N )
{
    double *data = (double *) calloc( N*N, sizeof(double) );
    double **M  = (double **) malloc( N  * sizeof(double*));

    for( int i = 0; i < N; i++ )
        M[i] = &data[i*N];

    return M;
}

/* Free's 2-D Contiguous Matrix */
void matrix_free( double ** M)
{
    free(M[0]);
    free(M);
}
// get time from both omp and mpi
double get_time(void)
{
#ifdef MPI
    return MPI_Wtime();
#endif

#ifdef OPENMP
    return omp_get_wtime();
#endif

    time_t time;
    time = clock();

    return (double) time / (double) CLOCKS_PER_SEC;
}
// error for command line parsing , need correct argument sor it wont run
void cli_error(void)
{
    printf("Please provide physical domain dimension as first argument to program\nand \'serial\' or \'parallel\' as the second argument, e.g.:\n\t$> ./cg 75 serial\n");
#ifdef MPI
    MPI_Finalize();
#endif
    exit(1);
}