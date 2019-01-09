#include "nbody_header.h"


void run_serial_problem(int nBodies, double dt, int nIters, char * fname, int threads)
{  // omp_set_num_threads(threads);
	// Open File and Write Header Info
	FILE * datafile = fopen(fname,"w");
	fprintf(datafile, "%+.*le %+.*le %+.*le\n", 10, (double)nBodies, 10, (double) nIters, 10, 0.0);

	// Allocate Bodies
	Body * bodies  = (Body *) calloc( nBodies, sizeof(Body) );
	randomizeBodies(bodies, nBodies, threads);

    double start = get_time();
   // for (int i = 0; i < nBodies; i++)
  //      fprintf(datafile, "%+.*le %+.*le %+.*le\n", 10, bodies[i].x, 10, bodies[i].y, 10, bodies[i].z);



	// Loop over timesteps
	for (int iter = 0; iter < nIters; iter++)
	{
		printf("iteration: %d\n", iter);

		// Output body positions to file
		for (int i = 0; i < nBodies; i++)
			fprintf(datafile, "%+.*le %+.*le %+.*le\n", 10, bodies[i].x, 10, bodies[i].y, 10, bodies[i].z);

		// Compute new forces & velocities for all particles
		compute_forces(bodies, dt, nBodies, threads);

		// Update positions of all particles
//#pragma omp parallel for num_threads(threads) private(i)
		for (int i = 0 ; i < nBodies; i++)
		{
			bodies[i].x += bodies[i].vx*dt;
			bodies[i].y += bodies[i].vy*dt;
			bodies[i].z += bodies[i].vz*dt;
		}

	}

	// Close data file
	fclose(datafile);

	double stop = get_time();

	double runtime = stop-start;
	double time_per_iter = runtime / nIters;
	long interactions = nBodies * nBodies;
	double interactions_per_sec = interactions / time_per_iter;

	printf("SIMULATION COMPLETE\n");
	printf("Runtime [s]:              %.3le\n", runtime);
	printf("Runtime per Timestep [s]: %.3le\n", time_per_iter);
	printf("Interactions per sec:     %.3le\n", interactions_per_sec);

	free(bodies);
}

// Randomizes all bodies to the following default criteria
// Locations (uniform random between -1.0 < r < 1.0 )
// Velocities (uniform random between -1.0e-3 < r < 1.0e3 )
// Masses (all equal at 1.0 / nBodies)
// You should make this more exotic
void randomizeBodies(Body * bodies, int nBodies, int threads)
{
	// velocity scaling term
	double vm = 1.0e-3;
//	omp_set_num_threads(threads);
//#pragma omp parallel for num_threads(threads) private(i) shared(vm,bodies,nBodies)
	for (int i = 0; i < nBodies; i++) {
        // Initialize position between -1.0 and 1.0
        if (i < nBodies / 2) {
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
        }  /*
             ///Initial
          bodies[i].x = 2.0 * (rand() / (double) RAND_MAX) - 1.0;
          bodies[i].y = 2.0 * (rand() / (double) RAND_MAX) - 1.0;
          bodies[i].z = 2.0 * (rand() / (double) RAND_MAX) - 1.0;
          // Intialize velocities
          bodies[i].vx = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
          bodies[i].vy = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
          bodies[i].vz = (2.0*vm * (rand() / (double)RAND_MAX) - vm);
*/
		// Initialize masses so that total mass of system is constant
		// regardless of how many bodies are simulated
		bodies[i].mass = 1.0 / nBodies;
	}
}

// Computes the forces between all bodies and updates
// their velocities accordingly
void compute_forces(Body * bodies, double dt, int nBodies, int threads)
{
	double G = 6.67259e-3;
	double softening = 1.0e-5;
  //  omp_set_num_threads(threads);
	// For each particle in the set
//#pragma omp parallel for num_threads(threads) private(i) shared(nBodies, dt, bodies, G, softening)
	for (int i = 0; i < nBodies; i++)
	{ 
		double Fx = 0.0;
		double Fy = 0.0;
		double Fz = 0.0;

		// Compute force from all other particles in the set
//#pragma omp parallel for num_threads(threads) private(j)
		for (int j = 0; j < nBodies; j++)
		{
			// F_ij = G * [ (m_i * m_j) / distance^3 ] * (location_j - location_i) 

			// First, compute the "location_j - location_i" values for each dimension
			double dx = bodies[j].x - bodies[i].x;
			double dy = bodies[j].y - bodies[i].y;
			double dz = bodies[j].z - bodies[i].z;

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
			double m_j = bodies[j].mass;
			double mGd = G * m_j / distance_cubed;
			Fx += mGd * dx;
			Fy += mGd * dy;
			Fz += mGd * dz;
		}

		// With the total forces on particle "i" known from this batch, we can then update its velocity
		// v = (F * t) / m_i
		// NOTE: as discussed above, we have left out m_1 from previous velocity computation,
		// so this is left out here as well
		bodies[i].vx += dt*Fx;
		bodies[i].vy += dt*Fy;
		bodies[i].vz += dt*Fz;
	}
}
