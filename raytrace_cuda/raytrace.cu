// cuda try implementation - Tomasz Michalik 05/22/18
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#define RADIUS 6.0
#define winMax 10

#define winY 10
#define PI 3.14159265359

typedef struct
{
    double x;
    double y;
    double z;
} vRay;


// 3D dot product
__device__
double dotproduct(vRay *V1, vRay *V2)
{
    return V1->x*V2->x + V1->y*V2->y + V1->z*V2->z;
}



// view
__device__
bool view_check(vRay * V, vRay * C, double * t)
{
    double VC = dotproduct(V, C);
    double CC = dotproduct(C, C);

    if ((VC*VC + RADIUS*RADIUS - CC) < 0)
        return false;

    *t = VC - sqrt(VC*VC +RADIUS*RADIUS - CC);

    return true;
}


// intersection of view and sphere
__device__
void intersect(vRay * I, vRay * V, double t)
{
    I->x = t*V->x;

    I->y = t*V->y;
    I->z = t*V->z;
}


// unit normal vector
__device__
void ray_normal(vRay *I, vRay *C, vRay *N)
{
    vRay ImC;

    ImC.x = I->x - C->x;
    ImC.y = I->y - C->y;
    ImC.z = I->z - C->z;

    double dot_ImC = dotproduct(&ImC, &ImC);

    N->x = (ImC.x) / sqrt(dot_ImC);
    N->y = (ImC.y) / sqrt(dot_ImC);
    N->z = (ImC.z) / sqrt(dot_ImC);
}


// create shadow ray and compute brightness
__device__
double brightness(vRay *I, vRay *L, vRay *N)
{
    vRay S;
    vRay LmI;

    LmI.x = L->x - I->x;
    LmI.y = L->y - I->y;
    LmI.z = L->z - I->z;

    double dot_LmI = dotproduct(&LmI, &LmI);

    S.x = (LmI.x) / sqrt(dot_LmI);
    S.y = (LmI.y) / sqrt(dot_LmI);
    S.z = (LmI.z) / sqrt(dot_LmI);

    // return the max between 0 and S.N
    double dots = dotproduct(&S, N);

    return (dots > 0) ? dots : 0.0;
}


// kernel function
__global__
void sampAlgo(double *grid, int n, int n_rays)
{
    // use cuda's random number generator
    int i =  blockDim.x*blockIdx.x + threadIdx.x;

    curandState_t state;
    curand_init(i, 0, 0, &state);

    double delta = (((double)n) / ((double)(2*winMax)));

    // set up light source and sphere center position
    vRay L;
    L.x = 4; L.y = 4; L.z = -1;

    vRay C;
    C.x = 0; C.y = 12; C.z = 0;

    vRay W; vRay V; vRay I; vRay N;

    for (int i = 0; i < n; ++i)
    {
        double t, theta, phi, b;

        do
        {
            phi = (double) curand_uniform(&state) * (double) M_PI;
            theta = (double) curand_uniform(&state)  * (double) M_PI;

            V.x  = sin(theta) * cos(phi);
            V.y  = sin(theta) * sin(phi);
            V.z  = cos(theta);

            W.x = (winY / V.y) * V.x;
            W.y = (winY / V.y) * V.y;
            W.z = (winY / V.y) * V.z;

        } while ((!view_check(&V, &C,  &t)) || (fabs(W.x) > winMax) || (fabs(W.z) > winMax));
// change to && , || produces similar timings
        intersect(&I, &V, t);

        ray_normal(&I, &C, &N);

        b = brightness(&I, &L, &N);

        double x = (W.x + (double)winMax);
        double z = (W.z + (double)winMax);
        x = x * delta;
        z = z * delta;

        grid[(int)x*n + (int)z] += b;

    }
}



int main(int argc, char **argv)
{
    // arg[1] = number of rays, arg[2] = number of grid points
    if (argc != 3)
    {
        printf("Invalid number of arguments.\n");
        exit(1);
    }

    int n = atoi(argv[1]);
    int n_rays = atoi(argv[2]);
    struct timeval start, end;
    gettimeofday(&start, NULL);

    srand(time(NULL));

    // allocate window (n x n)
    double * grid = (double *) calloc(n*n, sizeof(double));

    // copy data over launch kernel

    // Cuda malloc
    double * cuda_grid;
    cudaError_t _e;
    _e = cudaMalloc((void**)&cuda_grid, n*n * sizeof(double));
    if (_e != cudaSuccess)
        printf("Report Cuda error: %s\n", cudaGetErrorString(_e));

    //transfer  to gpu
    _e = cudaMemcpy(cuda_grid, grid, n*n*sizeof(double), cudaMemcpyHostToDevice);
    if (_e != cudaSuccess)
        printf("Report Cuda error: %s\n", cudaGetErrorString(_e));

    // run kernel
    int block_size =; // modify here
    int rays_per_thread = 10; // set default value
    int n_blocks = (n_rays + block_size - 1) / (block_size*rays_per_thread);
    printf("threads per block:%d number of blocks = %d\n", block_size, n_blocks);
    // blocks, threads per block parameters setx
    sampAlgo<<< n_blocks, block_size>>>(cuda_grid, n, rays_per_thread);
    _e = cudaGetLastError();

    //get from gpu
    _e = cudaMemcpy(grid, cuda_grid, n*n*sizeof(double), cudaMemcpyDeviceToHost);
    if (_e != cudaSuccess)
        printf("Cuda error: %s\n", cudaGetErrorString(_e));

    gettimeofday(&end, NULL);
    double m = 1000000;
    double t = ((end.tv_sec*m + end.tv_usec) - (start.tv_sec*m + start.tv_usec));
    printf("Cuda,%d, %g\n", n_rays, t / m);

    // write grid to file
    FILE * file = fopen("rays.out", "wb");

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            fwrite(&(grid[i*n + j]), sizeof(double), 1, file);
        }


    fclose(file);

    free(grid);

    return 0;
}