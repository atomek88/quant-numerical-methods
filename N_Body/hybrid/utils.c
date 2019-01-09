#include "nbody_header.h"
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

void print_inputs(long nBodies, double dt, int nIters, int nthreads )
{
	int mype = 0;
	int nprocs = 1;

	#ifdef MPI
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif

	if( mype == 0 )
	{
		printf("INPUT PARAMETERS:\n");
		printf("N Bodies =                     %ld\n", nBodies);
		printf("Timestep dt =                  %.3le\n", dt);
		printf("Number of Timesteps =          %d\n", nIters);
		printf("Number of Threads per Rank =   %d\n", nthreads);
		#ifdef MPI
		printf("Number of MPI Ranks =          %d\n", nprocs);
		#endif
		printf("BEGINNING N-BODY SIMLUATION\n");
	}
}
