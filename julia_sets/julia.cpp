//
// Created by Tomasz Michalik


#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include <vector>
#include <mpi.h>

#define MAX 1000
#define CX -0.7
#define MASTER 0
#define CY 0.26
int** allocMatrix(int rows, int cols);
int juliaCalc(double zx, double zy, double creal, double cimag);
void saveMatrix(int* T_0, int N);
void mpiPrint(int** Tmain, int row, int chunk, int length, MPI_File file, MPI_Status status);
int main(int argc, char* argv[]) {
    int i,j,k,c, N, nprocs, myrank, chunk, xstart, xend;
    double p, creal, cimag, dx, dy,zx, zy, start, end, time;
    MPI_Status status;
    MPI_File file;
    MPI_Datatype num_as_string;
    MPI_Datatype localarray;
    int**arr; int questions[2]; int coords[2]; int answer[2];
    int*pixl_line;
    double question;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    N = atoi(argv[1]);
    c = atoi(argv[2]); // arg for static = 0, dynamic = 1; *probs with io on midway linker
    creal = -0.7;
    cimag = 0.26;
    dx = 3.0/N; // x and y step resolution , physical space domain / N arr size
    dy = 2.0/N;
    //chunk = atoi(argv[3]); // chunk size parameter for dynamic loading only
    // test small array to pass
    if (N % nprocs != 0) {
        printf("Matrix not even, scheme to remainder\n");
    }
    pixl_line = (int*) malloc(N * sizeof(int));
    // master - workers
    MPI_Barrier(MPI_COMM_WORLD);
    if (c ==1){           // test with 4 proc runs
        if (myrank == MASTER) {
            start = MPI_Wtime(); // start time
            zx = -1.5;
            for (i = 0; i < N; i++) {
                // not synchronized, send, receive
                MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
                         &status); // answer[1] tells me what row it is, offset for printing or rmbr to save matrix data
                if (answer[1] >= 0) {  // will be N or ideal chunk size
                    //for (j = 1; j<= nprocs; j++) {
                    MPI_Recv(pixl_line, N, MPI_INT, answer[0], 2, MPI_COMM_WORLD, &status);
                    // arr[answer[0]][0]
                    saveMatrix(pixl_line, N);
                    //}
                }
                MPI_Send(&i, 1, MPI_INT, answer[0], 3, MPI_COMM_WORLD); // send info to each rank
                //printf("master sending %d th work to %d \n", i, answer[0]);
                MPI_Send(&zx, 1, MPI_DOUBLE, answer[0], 4, MPI_COMM_WORLD);
                zx = -1.5 + (dx * i);
            }
            // remaining terms ?
            k = -1; //end terms
            for (int i = 1; i < nprocs; ++i) {
                MPI_Recv(&answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                // puts in right spot -receiving here but not in order, may change to gather or sync to receive from answer[0] == 1 ... nprocs
                MPI_Recv(pixl_line, N, MPI_INT, answer[0], 2, MPI_COMM_WORLD, &status);
                //printf("master received from %d answer %d and %d \n", answer[0], arr[answer[1]][0], arr[answer[1]][1]);
                saveMatrix(pixl_line, N);
                // send termination signal
                MPI_Send(&k, 1, MPI_INT, answer[0], 3, MPI_COMM_WORLD);
                //printf("master sending termination signal to %d \n", answer[0]);
            }
            // end time
        }
        else {
            int i;
            answer[0] = myrank;
            answer[1] = -1; // for signal
            // send rank and message to master
            MPI_Send(answer, 2, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            while (true){ // loop for receiving work from master
                MPI_Recv(&i, 1, MPI_INT, MASTER, 3, MPI_COMM_WORLD, &status);
                if (i<0){ printf("My rank %d, received termination signal", myrank);
                    break;} // termination signal here
                answer[1] = i;

                MPI_Recv(&question, 1, MPI_DOUBLE, MASTER, 4, MPI_COMM_WORLD, &status);
                zy = -1.0;
                for (int j = 0; j<N; j++){
                    pixl_line[j] = juliaCalc(question, zy,creal, cimag);
                    zy = -1.0 + (dy *j);
                }
                MPI_Send(&answer, 2, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
                MPI_Send(pixl_line, N, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
                //printf("%d rank sending answer %d  %d tomasta \n", myrank, pixl_line[0], pixl_line[1]);
                // test send 1 int back greater than my rank
            }
        }
    } if (c == 0) {
        int remainder = N % nprocs;
        chunk = N / nprocs; // check if remainder - floor
        xstart = myrank * chunk;
        xend = chunk + xstart;
        if (remainder != 0 && myrank == (nprocs - 1)) {
            chunk = chunk + remainder;
            printf("last rank of %d has added remainder %d\n", myrank, remainder);
        }
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

            }
            // printf("end loop %d %d \n", i, j);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (i=0; i<nprocs;i++) {
            if (myrank == i) {/*
            printf("My rank %d data at startrow %d \n", myrank, xstart);
            for (i = 0; i < chunk; i++) {
                printf("\nrow:%d\n", i);
                for (j = 0; j < N; j++) {
                    printf("%d ", arr[i][j]);
                }
                printf("\n");
            }
        }
       */
                mpiPrint(arr, xstart, chunk, N, file, status);
            }
        }
    } else { printf("wrong parameters: choose (1) for dynamic, (0) for static");}
    MPI_Barrier(MPI_COMM_WORLD);
    // test results else {

    end = MPI_Wtime();
    printf("Time elapsed =%10.4f s\n", (end-start));
    // terminate MPI
    MPI_Finalize();

    return 0;
}

// Static printing** mpi print func for vector in order, chunk is # of rows in this print job
void mpiPrint(int** Tmain, int startrow, int chunk, int length, MPI_File file, MPI_Status status) {
    char *const fmt="%d ";
    char *const endfmt="%d\n";
    MPI_Datatype num_as_string;
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
    MPI_File_open(MPI_COMM_WORLD, "juliaDecomp.out",
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
    //  }
    return i;
}
// matrix allocation per proc
int** allocMatrix(int rows, int cols){
    int i;
    int* a = (int*) malloc(rows*cols *sizeof(int));
    int** aa = (int**) malloc(rows * cols * sizeof(int));
    for (i = 0; i< rows; i++){
        aa[i] = &(a[i*cols]); // each row pointer gets a col of data
    }
    return aa;
}
// save matrix
void saveMatrix(int* T_0, int N) {
    FILE *fp = fopen("mslave.txt", "a");
  //  file.open("/Users/tomaszmichalik/HPC/Prob3/mslave.txt", std::ios::out | std::ios_base::app); //std::ios::binary |

    for (int i = 0; i < N; i++) {
        int val = T_0[i];
        fprintf(fp, "%d ", val);
        //file << val << " ";
    }
    //file << "\n" << std::endl;
    fclose(fp);
    fprintf(fp, "\n");
    //file.close();
    //return;
}