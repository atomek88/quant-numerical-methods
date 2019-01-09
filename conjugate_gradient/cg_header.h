//
// Created by Tomasz Michalik on 5/17/18.
//

#ifndef PROB4_CG_HEADER_H
#define PROB4_CG_HEADER_H


#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <assert.h>
#include <vector>
//#ifdef MPI
#include <mpi.h>
//#endif
// main.c
void run_dense(long n);
void run_sparse(long n);
void run_parallel_sparse(long n, int mype, int nprocs, MPI_Status status);

// serial.c
void cg_dense(double ** A, double * x, double * b, long N);
void cg_sparse_poisson(double * x, double * b, long N);
void matvec( double * v, double ** M, double * w, long N );
void matvec_OTF( double * v, double * w, long N );
double dotp( double * a, double * b, long N);
void axpy( double alpha, double * w, double beta, double * v, long N);
void fill_A(double ** A, long N);
double find_b(long i, long j, long n);
void fill_b(double * b, long N);

// parallel.c
void parallel_cg_sparse_poisson(double * x, double * b, long N, int mype, int nprocs, int mypart, MPI_Status status);
void parallel_matvec_OTF( double * v, double * w, long N, int mype, int nprocs, MPI_Status status, double* nbuf, double* sbuf );
double parallel_dotp( double * a, double * b, long N, int mype, int nprocs, int iter);
void parallel_axpy( double alpha, double * w, double beta, double * v, long N, int mype, int nprocs, int iter);
void parallel_fill_b(double * b, long N, int mype, int nprocs);

// utils.c
void print_matrix(double ** A, long N );
void print_vector(double * x, long N );
void save_vector(double * x, long N, char * fname );
double ** matrix( long N );
void matrix_free( double ** M);
double get_time(void);
void cli_error(void);
void print_vector_p(double * x, long N, int rank, int nprocs);
void save_vector_p(double * x, long N, int rank, int nprocs, char * fname );
void save_vector_mpi(double * x, long N, int rank, int nprocs, MPI_File file, MPI_Status status);
#endif //PROB4_CG_HEADER_H
