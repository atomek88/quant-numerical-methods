//
// Created by Tomasz Michalik on 5/15/18.
//

#include <cstdio>
#include <cmath>
#include <ctime>
#include <string>
#include <ostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <complex>
using std::complex;
#define MAX 1000
#define CX -0.7
#define MASTER 0
#define CY 0.26

int** allocMatrix( int rows, int N);
int juliaCalc(double zx, double zy, double creal, double cimag);
int** setupArrIndex(int start, int ch, int Y);
void saveMatrix(int** T_0, int chunk, int N);
void mpiPrint(int** Tmain, int row, int chunk, int N, MPI_File f, MPI_Status stat);
int main(int argc, char* argv[]) {
    int i, j, k, c, N, nprocs, myrank, chunk, xstart, xend;
    double p, creal, cimag, dx, dy, zx, zy, start, end, time;
    MPI_Comm comm;
    MPI_Status status;
    int **arr;
    int questions[2];
    int answer[2];
    double question;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_File file;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    // std::ofstream myfile;
    N = atoi(argv[1]);
    c= atoi(argv[2]); // temp for y coord for testing
   // arg for static = 0, dynamic = 1; *probs with io on midway linker
    //chunk = atoi(argv[3]); // chunk size parameter for dynamic loading only
    printf("Rank:%d/%d starting with %dx%d matrix", myrank, nprocs, N, N);
    auto remainder = N % nprocs;
    if (remainder != 0) {
        printf("Matrix not evenly divisible among procs, scheme to remainder\n");
    }
    creal = -0.7;
    cimag = 0.26;
    dx = 3.0 / N; // x and y step resolution , physical space domain / N arr size
    dy = 2.0 / N;
    // divide matrix evenly with remainder for julia set calculation
    // try to use scatter to give each proc their coords and dim of their matrix to work on
    if (c == 0) {
        int remainder = 0;
        chunk = N / nprocs; // check if remainder - floor
        xstart = myrank * chunk;
        xend = chunk + xstart;
        if (remainder != 0 && myrank == (nprocs - 1)) {
            chunk = chunk + remainder;
            printf("last rank of %d has added remainder %d\n", myrank, remainder);
        }


        printf("rank %d : my submatrix coords are %d %d and chunksize is %d \n", myrank, xstart, xend, chunk);

        arr = allocMatrix(chunk, N);
        //arr = setupArrIndex(xstart, chunk, N);
        // try to divide up all work, then print from 0 - - n into file
        MPI_Barrier(MPI_COMM_WORLD);
        // run function for my subsection of array

        start = MPI_Wtime();
        for (i = 0; i < chunk; i++) { //xstart; i < xend;++i

            zx = -1.5 + (dx * (i + xstart));
            for (j = 0; j < N; j++) {

                zy = -1.0 + (dy * j);
                //printf("rank:%d here check initial value %d \n",myrank,  arr[i][j]);
                k = juliaCalc(zx, zy, creal, cimag);
                arr[i][j] = k;
                // printf("%d ", k);
                //printf("rank:%d here my values %d \n",myrank,  arr[i][j]);

            }
            // printf("end loop %d %d \n", i, j);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        /*
        for (i=0; i<nprocs;i++) {
        if (myrank == i) {
          //  printf("My rank %d data at startrow %d \n", myrank, xstart);
            for (i = 0; i < chunk; i++) {
                printf("\nrow:%d\n", i);
                for (j = 0; j < N; j++) {
                    printf("%d ", arr[i][j]);
                }
                printf("\n");
            }
        } */
        MPI_Barrier(MPI_COMM_WORLD);
        mpiPrint(arr, xstart, chunk, N, file, status);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        // test results



        MPI_Finalize();
        printf("Time elapsed =%10.4f s\n", (end - start));
        return 0;
    }



// julia
int juliaCalc(double zx, double zy, double creal, double cimag){
    double tmp; int i;
    double result;
    double zreal; double zimag;
    zreal = zx; zimag = zy;
    result = zreal * zreal + zimag * zimag;
    // try this implementation : other () double ax = 0.0; double ay = 0.0;
    //while ( (iteration <= MAX) && (result > 4.0) ) {
    for (i = 0; i <= MAX; i++) {
        tmp = (zreal * zreal) - (zimag * zimag);
        zimag = 2 * zreal * zimag + cimag;
        zreal = tmp + creal;
        result = (zreal * zreal) + (zimag * zimag);

        if (result > 4.0) break;
    }
    //}
    return i;
}
//
// mpi print func for vector in order, chunk is # of rows in this print job
void mpiPrint(int** Tmain, int startrow, int chunk, int length, MPI_File file, MPI_Status status) {
    char *const fmt="%d ";
    char *const endfmt="%d\n";
    MPI_Datatype num_as_string;
    printf("start row %d, chunk %d, length %d printing\n", startrow, chunk, length);
    MPI_Datatype localarray;
    const int charspernum=4;
    MPI_Type_contiguous(charspernum, MPI_CHAR, &num_as_string);
    MPI_Type_commit(&num_as_string);
    char *data_as_txt = (char*) malloc(chunk * length * charspernum*sizeof(char));
    int count = 0;
    for (int i=0; i<chunk; i++) {
        for (int j=0; j<length-1; j++) {
            sprintf(&data_as_txt[count*charspernum], fmt, Tmain[i][j]);
            count++;
        }
        sprintf(&data_as_txt[count*charspernum], endfmt, Tmain[i][length-1]);
        count++;
    }
    int globalsizes[2] = {length, length};
    int localsizes [2] = {chunk, length};
    int starts[2]      = {startrow, 0};
    MPI_Type_create_subarray(2, globalsizes, localsizes, starts,
                             MPI_ORDER_C, num_as_string,
                             &localarray);
    MPI_Type_commit(&localarray);

    MPI_File_open(MPI_COMM_WORLD, "mpijFull.txt",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);
    //MPI_Offset displace = rank * length * sizeof(int); // for static decomp can add row
    MPI_File_set_view(file, 0, MPI_INT, localarray,
                      "native", MPI_INFO_NULL);
/* Change here to print 'Tmain' to print Binary. All above related to char conversion can be commented out */
    MPI_File_write_all(file, Tmain, length*chunk, num_as_string, &status); // diff between write all and write?
    MPI_File_close(&file);
    MPI_Type_free(&localarray);
    MPI_Type_free(&num_as_string);

}
int** setupArrIndex(int start, int chunk, int Y){
    printf("rank: setting up ranks index subarray, size:%d \n", chunk );
    int i,j;
    int end = start + chunk;
    int* a = (int*) malloc(Y *sizeof(int));
    int** aa = (int**) malloc(chunk*Y * sizeof(int));
    for (i=start; i<=end; i++){
        aa[i] = &(a[i*Y]);}
    return aa;
}
// matrix allocation per proc
int** allocMatrix(int rows, int N){
    int i;
    int* a = (int*) malloc(N *rows*sizeof(int)); // need to add rows * N?
    int** aa = (int**) malloc(rows * N * sizeof(int));
    for (i = 0; i< rows; i++){
        aa[i] = &(a[i*N]); // each row pointer gets a col of data
    }
    return aa;
}
void saveMatrix(int** T_0, int chunk, int N) {
    printf("Trying to print something.. %d , %d \n", T_0[0][0], T_0[1][0]);
    std::ofstream file;
    file.open("/Users/tomaszmichalik/HPC/Prob3/outdecompT.txt", std::ios::out | std::ios_base::app); //std::ios::binary |

    for (int i = 0; i < chunk; i++) {
        for (int j = 0; j < N; j++) {
            int val = T_0[i][j];
            printf("%d ",val);
            file << val << " ";
        }//ofile.write((char*) &f, sizeof(float));
        file << "\n" << std::endl;
        printf("\n");
    }
    file.close();
    //return;
}