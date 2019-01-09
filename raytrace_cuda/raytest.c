//
// Created by Tomasz Michalik on 6/1/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#define RADIUS 6.0
#define winMax 10
#define PI 3.14159265359
#define winY 10


// free grid
void freegrid(double **arr)
{
    free(arr[0]);
    free(arr);
}


typedef struct
{
    double x;
    double y;
    double z;
} vRay;


// check if unit vector
bool unit_vector(vRay *V)
{
    int norm = sqrt(V->x*V->x + V->y*V->y + V->z*V->z);

    if (norm == 1) return true; else return false;
}


// 3D dot product
double dotproduct(vRay *V1, vRay *V2)
{
    return V1->x*V2->x + V1->y*V2->y + V1->z*V2->z;
}

// window intersection scalar
// assumes V is unit vector and returns true if scalar exists
bool view_check(vRay * V, vRay * C, double * t)
{
    double VC = dotproduct(V, C);
    double CC = dotproduct(C, C);

    if ((VC*VC + RADIUS*RADIUS - CC) < 0)
        return false;

    *t = VC - sqrt(VC*VC + RADIUS*RADIUS - CC);

    return true;
}

// intersection of ray and sphere
void intersect(vRay * I, vRay * V, double t)
{
    I->x = t*V->x;
    I->y = t*V->y;
    I->z = t*V->z;
}


// unit normal vector of sphere
void normalize_sphere(vRay *I, vRay *C, vRay *N)
{
    vRay * Im = malloc(sizeof(vRay));

    Im->x = Im->x - C->x;
    Im->y = Im->y - C->y;
    Im->z = Im->z - C->z;

    double dot_I = dotproduct(Im, Im);

    N->x = (Im->x) / sqrt(dot_I);
    N->y = (Im->y) / sqrt(dot_I);
    N->z = (Im->z) / sqrt(dot_I);

    free(Im);
}


// create shadow ray and compute brightness
double brightness(vRay *I, vRay *L, vRay *N)
{
    vRay * S = malloc(sizeof(vRay));
    vRay * LmI = malloc(sizeof(vRay));
    LmI->x = L->x - I->x;
    LmI->y = L->y - I->y;
    LmI->z = L->z - I->z;

    double dotL = dotproduct(LmI, LmI);

    S->x = (LmI->x) / sqrt(dotL);
    S->y = (LmI->y) / sqrt(dotL);
    S->z = (LmI->z) / sqrt(dotL);

    // return the max between 0 and S.N
    double dots = dotproduct(S, N);

    free(S);
    free(LmI);
    if (dots > 0) { return dots;}
    return 0.0;
}
// matrix allocation per proc
double** allocMatrix(int rows, int cols){
    int i;
    double* a = (double*) calloc(rows*cols ,sizeof(double));
    double** aa = (double**) malloc(rows * sizeof(double*));
    for (i = 0; i< rows; i++){
        aa[i] = &(a[i*cols]); // each row pointer gets a col of data
    }
    printf("Allocating matrix with 0s\n");
    for (i=0; i< rows; i++){
        for( int j=0; j<cols; j++) {
            aa[i][j] = 0.0;
        }
    }
    return aa;
}

// float((double)rand()/(double)(RAND_MAX/a))
double rand_d(double upper_limit)
{
    return ((double)rand()/(double)(RAND_MAX/upper_limit));
}
// randnum
double randNum(int min, int max){
   // printf("Rand numb genreator: %f\n", r);
    return ((double)rand()/(double)RAND_MAX * (max - min) + min);
}

int main(int argc, char **argv)
{
    struct timeval start, end;
    gettimeofday(&start, NULL);


    // arg[1] = number of rays, arg[2] = number of grid points
    if (argc != 3)
    {
        printf("Invalid number of arguments.\n");
        exit(1);
    }
    int n = atoi(argv[1]);
    int rays = atoi(argv[2]);


    srand(time(NULL));

    // allocate window (nxn)
    double ** grid = allocMatrix(n,n);

    double delta = (((double)n) / ((double)(2*winMax)));

    // set L
    vRay * L = malloc(sizeof(vRay));
    L->x = 4; L->y = 4; L->z = -1;

    vRay * C = malloc(sizeof(vRay));
    C->x = 0; C->y = 12; C->z = 0;

    vRay * W = malloc(sizeof(vRay));
    vRay * V = malloc(sizeof(vRay));
    vRay * I = malloc(sizeof(vRay));
    vRay * N = malloc(sizeof(vRay));

    for (int i = 0; i < rays; ++i)
    {
        double t, theta, phi, b;

        do
        {
            theta = randNum(0,PI); // or change to just PI - PI
            phi   = randNum(0,PI);
           // double temp = sqrt(1 - pow(cos(theta), 2)); // trig identities
           // V->x = temp * cos(phi);
            //V->y = temp * sin(phi);
            //V->z = cos(theta);
            V->x  = sin(theta) * cos(phi);             // known works
            V->y  = sin(theta) * sin(phi);
            V->z  = cos(theta);

            W->x = (winY / V->y) * V->x;
            W->y = (winY / V->y) * V->y;
            W->z = (winY / V->y) * V->z;

        } while ((!view_check(V, C,  &t)) || (fabs(W->x) > winMax) || (fabs(W->z) > winMax));

        intersect(I, V, t);

        normalize_sphere(I, C, N);

        b = brightness(I, L, N);

        double x = (W->x + (double)winMax);
        double z = (W->z + (double)winMax);
        x = x * delta;
        z = z * delta;

        grid[(int)x][(int)z] += b;
    }

    // print execution time
    gettimeofday(&end, NULL);
    double m = 1000000;
    double t = ((end.tv_sec*m + end.tv_usec) - (start.tv_sec*m + start.tv_usec));
    printf("serial,%d, %g\n", rays, t / m);


    // write grid to file
    FILE * file = fopen("rays.out", "wb");

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            fwrite(&(grid[i][j]), sizeof(double), 1, file);

    fclose(file);

    freegrid(grid);
    free(I);
    free(N);
    free(W);
    free(L);
    free(C);

    return 0;
}