// ping transfer communication timing between nodes
// Tomasz MIchalik ping pong bandwidth mpi test
// start with 250 ints for kb, then * 1000 for each iteration

// INCOMPLETE, functions decompose into main to make work - see pingerA
#include <cstdio>
#include <cmath>
#include <ctime>
#include <string>
//#include <ostream>
#include <cstdlib>

#include <vector>
#include <mpi.h>

#define A 0
#define B 1
#define ping 10
#define pong 10
#define TIMES 100 // increase
#define LEN 250000000 // add 1 for overflow prot?

void processA(int* s, int* r, int size);
void processB(int* s, int* r, int size);

int main(int argc, char* argv[]) {
    int i, j, nprocs, myrank, N;
    int* sendmsg; int* rcvmsg;
    MPI_Status status;
    double bandw;
   // N = atoi(argv[1]); // N set for KB, MB, GB
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    // check proc number
   // MPI_Barrier(MPI_COMM_WORLD);
    int size = 250;
    for (int y = 1; y<4; y++) {
        if (myrank == A) {
            //printf("A processing\n");
            processA(sendmsg, rcvmsg,  size);
        } else if (myrank == B) {
            //printf("B processing\n");
            processB(sendmsg, rcvmsg, size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
// remove procs, all same
// process A
void processA(int* send, int* rec){
    int i, N, size;
    MPI_Status status;
    double start, finish, time, bandw; // timing for sending messages
    size = 250; // start with kb 250 ints = ~ 1000 bytes (not 1024)
    send = (int*) malloc(LEN*sizeof(int));
    rec = (int*) malloc(LEN*sizeof(int));
    N = sizeof(int);


    //MPI_Barrier(MPI_COMM_WORLD);
    // set message size
    // 3 time iteration for kb, mb, gb
    for (i= 1; i< 4; i++) {
        time = 0.00000000000;
        printf("\nRound trip test for: \n");
        printf("Message length in # of ints= %d\n", size);
        printf("Message size = %d bytes\n", N * size);
        printf("Num reps = %d\n", TIMES);
        start = MPI_Wtime();
        // loop for sending messages, try 1000 first
        for (i = 1; i <= TIMES; i++) {
            // use mpi ssend to synchronously await the receive end
            MPI_Ssend(send, size, MPI_INT, B, ping, MPI_COMM_WORLD);
            // ping
            MPI_Recv(rec, size, MPI_INT, B, pong, MPI_COMM_WORLD, &status);
        } // save message passing?
        finish = MPI_Wtime();
        time = finish - start;
        bandw = 2.0 * TIMES * N * size / time; // bandwidth is 2 ranks * number * size / time
        printf("Bandwidth %d= %f Kyte/sec\n", size, bandw * 8.0/ 1000.0);
        printf("          = %f Mbit/sec\n", bandw * 8.0 / 1000000.0);
        printf("          = %f Gbit/sec\n", bandw * 8.0 / 1000000000.0);
    //    printf("%d\t  %f    %f\n", size, time / 200., (double) (2 * snum * 100 * size) / time);
        size = 1000*size;
    }

}
// process B
void processB(int* send, int* rec, int siz){
    int i;
    MPI_Status status;
    send = (int*) malloc(LEN*sizeof(int));
    rec = (int*) malloc(LEN*sizeof(int));
    int size =siz;
// 3 times for kb, mb, gb

    for (i =1; i<= TIMES; i++){

            MPI_Recv(rec, size, MPI_INT, A, ping, MPI_COMM_WORLD, &status);
            MPI_Ssend(send, size, MPI_INT, A, pong, MPI_COMM_WORLD);
             }
    printf("iteration %d size b: %d",i, size);

}
/* barrier example
 * MPI_Comm_size(MPI_COMM_WORLD, &size);
for ( int ii = 0; ii < size; ++ii ) {
    if ( rank == ii ) {
        // my turn to write to the file
        writeStuffToTheFile();
    }
    MPI_Barrier(MPI_COMM_WORLD);*/