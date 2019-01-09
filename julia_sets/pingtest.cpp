// ping transfer communication timing between nodes
// Tomasz Michalik MPCS 51087
#include <cstdio>
#include <cmath>
#include <ctime>

#include <cstdlib>

#include <mpi.h>

#define A 0
#define B 1
#define TIMES 10000 //10000

#define ping 10
#define pong 10

int main(int argc, char* argv[]) {
    int i, j, nprocs, myrank, N, length, ierr;
    int* sendmsg; int* rcvmsg; double* arr;
    double start, end, time, bandw;
    MPI_Status status;
    N = sizeof(int); // set to size of int

    ierr = MPI_Init(&argc, &argv);
    ierr|=MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ierr|=MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    // check proc number
    if (nprocs !=2){ printf("Only valid with 2 mpi ranks");
        exit(0);}

    length = 1;
    sendmsg = (int*) malloc(sizeof(int));
    rcvmsg = (int*) malloc(sizeof(int));
    //sendmsg = (int*) malloc(LEN*sizeof(int));
   //rcvmsg = (int*) malloc(LEN*sizeof(int));
    // barrier
    ierr = MPI_Barrier(MPI_COMM_WORLD);

    arr = (double*) malloc(4*sizeof(int));
    if (myrank == A) {
        for (i=1; i<=4;  i++) {
        int r = 6;
        sendmsg = &r;
            time = 0.000000000;
            // Roundtrip testing timing
            printf("\nRound trip test for: \n");
            printf("Message length = %d\n", length);
            printf("Message size = %d bytes\n", N * length);
            printf("Num reps = %d\n", TIMES);
            // start times
            start = MPI_Wtime();
            for (int j = 1; j <= TIMES; j++) {  // use Ssend?
                ierr = MPI_Ssend(sendmsg, length, MPI_INT, B, ping, MPI_COMM_WORLD);
                ierr = MPI_Recv(rcvmsg, length, MPI_INT, B, pong, MPI_COMM_WORLD, &status);
            }
            end = MPI_Wtime();
            // get average times
            time = end - start;
            arr[i-1] = time/TIMES;
            printf("Round Trip avg = %.2lf uSec\n", time / TIMES);
            /*
            bandw = 2.0 * TIMES * N * length / time; // bandwidth is 2 ranks * number * size / time
            printf("Bandwidth = %f Byte/sec\n", bandw);
            printf("          = %f Mbit/sec\n", bandw * 8.0 / 1000000.0);
            printf("          = %f Gbit/sec\n", bandw * 8.0 / 1000000000.0);

            length = 100 * length;  */
        }
        double avg;
        for (i =0; i< 4; i++){
            avg += arr[i];}

        printf("My average of averages: %f", avg/4);
    }
    else  { //if (myrank == B)
    int t = 2;
    sendmsg = &t;
        for (i=1; i<=4; i++){
            for (j = 1; j<= TIMES; j++) {
                ierr = MPI_Recv(rcvmsg, length, MPI_INT, A, pong, MPI_COMM_WORLD, &status);
                ierr = MPI_Ssend(sendmsg, length, MPI_INT, A, ping, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Finalize();
    return 0;
}
