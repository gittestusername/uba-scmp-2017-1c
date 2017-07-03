/*
  "Hello World" MPI Test Program
 */
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <chrono> // for high_resolution_clock.
#include <math.h>
#include <sstream>


#define BUFSIZE 128
#define TAG 0

using namespace std;

long double diff(long double[] &A, long double[] &B, int nX, int nY) {
    long double sum = 0;
    for (int i = 0; i <nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            sum += pow(B[i][j] - A[i][j], 2);
        }
    }
    return sqrt(sum);
}

void printMat(long double[] m) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m[0].size(); ++j) {
            cout << m[i][j] << " " ;
        }
        cout << endl;
    }
    cout << endl;
}



int main(int argc, char *argv[]) {
    char idstr[32];
    char buff[BUFSIZE];
    int numprocs;
    int myid;
    int i;
    MPI_Status stat;
    /* MPI programs start with MPI_Init; all 'N' processes exist thereafter */
    MPI_Init(&argc, &argv);
    /* find out how big the SPMD world is */
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    /* and this processes' rank is */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    /* At this point, all programs are running equivalently, the rank
       distinguishes the roles of the programs in the SPMD model, with
       rank 0 often used specially... */
    long double xMax = 2.0;
    long double yMax = 2.0;
    long double tMax = 0.05;
    long double nu = 0.1; //viscosidad
    long double rho = 1.0;  //densidad
    long double dx = 1.0 / 20.0;
    long double dy = 1.0 / 20.0;
    long double dt = 0.00001;
    int nX = round(xMax / dx) + 1;
    int nY = round(yMax / dy) + 1;
    int nT = round(tMax / dt) + 1;
    long double al = 0.5;
    bool upwind = false;
    long double fixedPointError = 0.000001;
    long double minFixedPointIters = 10;
    bool debug = true;
    long double xc = xMax / 2;
    long double yc = yMax / 2;
    long double rMax = 0.8 * min(xMax, yMax) / 2.0;
    long double rMin = 0.3 * min(xMax, yMax) / 2.0;
    long double fanTurns = 2.0;
    long double pi = atan(1) * 4;



    long double U0[BUFSIZE];
    long double U1[BUFSIZE];
    long double U2[BUFSIZE];
    long double V0[BUFSIZE];
    long double V1[BUFSIZE];
    long double V2[BUFSIZE];
    long double P0[BUFSIZE];
    long double P1[BUFSIZE];
    long double P2[BUFSIZE];



    if (myid == 0) {
        cerr << myid << "We have " << numprocs << " processors\n" << endl;
        for (i = 1; i < numprocs; i++) {
            sprintf(buff, "Hello %d! ", i);
            MPI_Send(buff, BUFSIZE, MPI_CHAR, i, TAG, MPI_COMM_WORLD);
        }
        for (i = 1; i < numprocs; i++) {
            MPI_Recv(buff, BUFSIZE, MPI_CHAR, i, TAG, MPI_COMM_WORLD, &stat);
            printf("%d: %s\n", myid, buff);
        }
    } else {
        /* receive from rank 0: */
        MPI_Recv(buff, BUFSIZE, MPI_CHAR, 0, TAG, MPI_COMM_WORLD, &stat);
        sprintf(idstr, "Processor %d ", myid);
        strncat(buff, idstr, BUFSIZE - 1);
        strncat(buff, "reporting for duty\n", BUFSIZE - 1);
        /* send to rank 0: */
        MPI_Send(buff, BUFSIZE, MPI_CHAR, 0, TAG, MPI_COMM_WORLD);
    }

    /* MPI programs end with MPI Finalize; this is a weak synchronization point */
    MPI_Finalize();
    return 0;
}