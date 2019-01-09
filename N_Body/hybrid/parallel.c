#include "nbody_header.h"

#ifdef MPI
#include "nbody_header.h"


void run_parallel_problem(int nBodies, double dt, int nIters, char * fname, int threads)
{  // Open File and Write Header Info
    int nprocs, myrank,mypart,rightn, leftn;
    MPI_Status status;
    MPI_File file;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    mypart = nBodies / nprocs;
    // datatyp[e = type
    MPI_Datatype type;
    MPI_Type_contiguous(7, MPI_DOUBLE, &type);
    MPI_Type_commit(&type);

    printf("Running parallel problem, rank:%d size:%d \n", myrank, mypart);
    omp_set_num_threads(threads);
	// Allocate Bodies
	Body* bodies;
    MPI_Request request[2];
    // distribute particles, randomize from master, other ranks receive
    if (myrank == 0) {
        printf("Master dividing bodies here\n");
        initialize_IO(nBodies, nIters, fname);
        for (int i = 1; i<nprocs; i++) {
            bodies  = (Body *) calloc( mypart, sizeof(Body) );
            printf("Master creating for rank%d with tag:%d \n", i, i+1);
            randomizeBodies(bodies, mypart, threads);
            MPI_Send(bodies, mypart, type, i, i+1, MPI_COMM_WORLD);
        } // master do its own work
        bodies  = (Body *) calloc( mypart, sizeof(Body) );
        randomizeBodies(bodies, mypart, threads);
        printf("Master creates its own particles! test:%f \n",bodies[0].x);

    } else {
        bodies  = (Body *) calloc( mypart, sizeof(Body) );
        printf("rank%d receiving bodies, tag:%d \n",myrank, myrank+1);
        MPI_Recv(bodies, mypart, type, 0, myrank+1, MPI_COMM_WORLD, &status);
        printf("rank%d received, last test:%f \n",myrank, bodies[0].x);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double start = get_time();
    Body * sendbuf  = (Body*) calloc( mypart, sizeof(Body) );
    Body * recbuf = (Body *) calloc(mypart, sizeof(Body) );

    // Loop over timesteps
    for (int iter = 0; iter < nIters; iter++) {
        rightn = myrank+1 >= nprocs ? 0 : myrank+1;
        leftn = myrank-1 < 0 ? nprocs -1 : myrank -1;
        printf("iteration: %d\n", iter);
        memcpy(sendbuf, bodies, mypart * sizeof(Body));

        // copy into send buf, write to file  **** had to hardcore nIters, for some reason changing to 0 after 1st iteration
        distributed_write_timestep(bodies, mypart, iter, 10, nprocs, myrank, file, status);
        // neighbors

        // compute on my own part
        compute_forces_multi_set(bodies, sendbuf, dt, mypart, threads);
        // interchange information loop, compute on everyons part
        for (int k =0; k < nprocs-1; k++) {

            printf("My rank:%d k=%d sending %d, receiving from:%d \n", myrank, k,leftn, rightn);
            MPI_Irecv(recbuf, mypart, type, rightn, 0, MPI_COMM_WORLD, &request[0]);
            MPI_Isend(sendbuf, mypart, type, leftn, 0, MPI_COMM_WORLD, &request[1]);
// do stuff, compute forces here
            MPI_Waitall(2, request, &status);
            MPI_Barrier(MPI_COMM_WORLD);
            printf("rank%d received from %d , send test:%f k=%d\n",myrank, rightn, recbuf[0].x, k);
            rightn = rightn + 1 >= nprocs ? 0 : rightn + 1;
            leftn = leftn - 1 < 0 ? nprocs - 1 : leftn - 1;
            MPI_Barrier(MPI_COMM_WORLD);
            printf("My rank:%d k=%d AFTER:: next send %d, next recv from:%d \n", myrank, k,leftn, rightn);
            compute_forces_multi_set(bodies, recbuf, dt, mypart, threads);
            recbuf = (Body *) calloc(mypart,sizeof(Body));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //compu
        // Update positions of all particles
#pragma omp parallel for num_threads(threads) private(i)
        for (int i = 0 ; i < mypart; i++)
        {
            bodies[i].x += bodies[i].vx*dt;
            bodies[i].y += bodies[i].vy*dt;
            bodies[i].z += bodies[i].vz*dt;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }


/*
    for (int i = 0; i < 2; i++) {
        printf("see if i can print over file, rank%d, i-%d\n", myrank, i);
        distributed_write_timestep(bodies, mypart, i, 2, nprocs, myrank, file, status);
    } */


    MPI_Barrier(MPI_COMM_WORLD);
    double stop = get_time();
    // timing
	double runtime = stop-start;
	double time_per_iter = runtime / nIters;
	long interactions = nBodies * nBodies;
	double interactions_per_sec = interactions / time_per_iter;
if (myrank == 0) {
	printf("SIMULATION COMPLETE\n");
	printf("Runtime [s]:              %.3le\n", runtime);
	printf("Runtime per Timestep [s]: %.3le\n", time_per_iter);
	printf("Interactions per sec:     %.3le\n", interactions_per_sec);
}
	free(bodies);
}

void compute_forces_multi_set(Body * local, Body * remote, double dt, int n, int threads)
{   double G = 6.67259e-3;
	double softening = 1.0e-5;
    printf("Compute forces loop\n");
	// For each particle in the set
	omp_set_num_threads(threads);
#pragma omp parallel for num_threads(threads) private(i) shared(nBodies, dt, bodies, G, softening)
	for (int i = 0; i < n; i++)
	{
		double Fx = 0.0;
		double Fy = 0.0;
		double Fz = 0.0;
#pragma omp parallel for num_threads(threads) private(j, remote, local)
		// Compute force from all other particles in the set
		for (int j = 0; j < n; j++)
		{
			// F_ij = G * [ (m_i * m_j) / distance^3 ] * (location_j - location_i)

			// First, compute the "location_j - location_i" values for each dimension
			double dx = remote[j].x - local[i].x;
			double dy = remote[j].y - local[i].y;
			double dz = remote[j].z - local[i].z;

			// Then, compute the distance^3 value
			// We will also include a "softening" term to prevent near infinite forces
			// for particles that come very close to each other (helps with stability)

			// distance = sqrt( dx^2 + dy^2 + dz^2 )
			double distance = sqrt(dx*dx + dy*dy + dz*dz + softening);
			double distance_cubed = distance * distance * distance;

			// Now compute G * m_2 * 1/distance^3 term, as we will be using this
			// term once for each dimension
			// NOTE: we do not include m_1 here, as when we compute the change in velocity
			// of particle 1 later, we would be dividing this out again, so just leave it out
			double m_j = local[j].mass;
			double mGd = G * m_j / distance_cubed;
			Fx += mGd * dx;
			Fy += mGd * dy;
			Fz += mGd * dz;
		}

		// With the total forces on particle "i" known from this batch, we can then update its velocity
		// v = (F * t) / m_i
		// NOTE: as discussed above, we have left out m_1 from previous velocity computation,
		// so this is left out here as well
		local[i].vx += dt*Fx;
		local[i].vy += dt*Fy;
		local[i].vz += dt*Fz;
	}
}

void parallel_randomizeBodies(Body * bodies, int nBodies_per_rank, int threads)
{// velocity scaling term
	double vm = 1.0e-3;
    printf("In parallel randomizing bodies\n");

	for (int i = 0; i < nBodies_per_rank; i++) {
  // Initialize position between -1.0 and 1.0
        if (i < nBodies_per_rank / 2) {
            bodies[i].x = 2.0 * (rand() / (double) RAND_MAX) - 1.0;
            bodies[i].y = 1.0 * (rand() / (double) RAND_MAX) - 1.0;
            bodies[i].z = -1.3 *(rand() / (double) RAND_MAX) - 1.0;
        } else {
            bodies[i].x = 2.3 * (rand() / (double) RAND_MAX) - 1.0;
            bodies[i].y = 1.4 * (rand() / (double) RAND_MAX) - 1.0;
            bodies[i].z = 2.0 * (rand() / (double) RAND_MAX) - 1.0;
        } // velocities change up
        if (bodies[i].x < 0){
            // velocities
            bodies[i].vx = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vy = (1.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vz = (-2.0*vm * (rand() / (double)RAND_MAX) - vm);
        }
        else if (bodies[i].x > 0){
            // velocities
            bodies[i].vx = (-2.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vy = (-1.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vz = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
        } else if (bodies[i].y > -0.3 && bodies[i].y < 0.3){
            // velocities
            bodies[i].vx = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vy = (-1.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vz = vm * (rand() / (double)RAND_MAX) - vm;
        } else {
            bodies[i].vx = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vy = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
            bodies[i].vz = (1.0*vm * (rand() / (double)RAND_MAX) - vm);
        }

		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated
		bodies[i].mass = 1.0 / nBodies_per_rank;
	}

}

// Opens MPI file, and writes header information (nBodies, iterations)
void distributed_write_timestep(Body * local_bodies, int nBodies_per_rank, int timestep, int nIters, int nprocs, int mype, MPI_File  fh, MPI_Status stat)
{  // print function in submatrix and convert ASCII format
    //    fprintf(datafile, "%+.*le %+.*le %+.*le\n", 10, bodies[i].x, 10, bodies[i].y, 10, bodies[i].z);
    char *const fmt="%.9le ";
    char *const endfmt="%.9le\n";
    MPI_Datatype num_as_string;
    int startrow = mype *nBodies_per_rank;
    MPI_Datatype localarray;
    const int charspernum=54;
    MPI_Type_contiguous(charspernum, MPI_CHAR, &num_as_string);
    MPI_Type_commit(&num_as_string); // create a one liner for each body of file, each particle
    char *data_as_txt = (char*) malloc(nBodies_per_rank *charspernum*sizeof(char));
    int count = 0;
    for (int i=0; i<nBodies_per_rank; i++) {//sprintf (buffer, "%d plus %d is %d", a, b, a+b);
        //    sprintf(&data_as_txt[count * charspernum], fmt, local_bodies[i].x, fmt, local_bodies[i].y, endfmt, local_bodies[i].z);
            sprintf(&data_as_txt[i * charspernum], "%+.*le %+.*le %+.*le\n",10, local_bodies[i].x, 10, local_bodies[i].y, 10, local_bodies[i].z);
           // count++;
    }
    int step = (nBodies_per_rank*nprocs) * timestep;
  // testprint  printf("%d: %s\n", mype, data_as_txt);
    int globalsizes[2] = {(int)nBodies_per_rank*nprocs*nIters, 0}; // * iters too
    int localsizes [2] = {(int)nBodies_per_rank, 0};
    int starts[2]      = {step+startrow, 0}; // need row to start across each timestep, stride it
    printf("Global: %d, rank%d, startrow:%d, step:%d local %d iter:%d\n", globalsizes[0], mype,startrow, step,localsizes[0], timestep);
    MPI_Type_create_subarray(1, globalsizes, localsizes, starts,
                             MPI_ORDER_C, num_as_string,
                             &localarray);
    MPI_Type_commit(&localarray);

    MPI_File_open(MPI_COMM_WORLD, "parCorr.txt",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_CHAR, localarray,
                      "native", MPI_INFO_NULL);
/* Change here to print 'Tmain' to print Binary. All above related to char conversion can be commented out */
    MPI_File_write_all(fh, data_as_txt, (int)nBodies_per_rank, num_as_string, &stat); // diff between write all and write?
    MPI_File_close(&fh);
    MPI_Type_free(&localarray);
    MPI_Type_free(&num_as_string);
}

void initialize_IO( long nBodies, long nIters, char* fname){
    FILE * datafile = fopen(fname,"w");
    fprintf(datafile, "%+.*le %+.*le %+.*le\n", 10, (double)nBodies, 10, (double) nIters, 10, 0.0);
    /// void mpiPrint(int** Tmain, int startrow, int chunk, int length,  MPI_Status status)

}
#endif