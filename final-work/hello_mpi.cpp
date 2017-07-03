#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <chrono>       // for high_resolution_clock.
#include <math.h>
#include <vector>
#include <utility>      // std::pair
#include <sstream>
#include "mat2.h"
#include "parameters.h"


#define TAG 0

using namespace std;


long double diff(mat2 A, mat2 B) {
    long double sum = 0;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            sum += pow(B.at(i, j) - A.at(i, j), 2);
        }
    }
    return sqrt(sum);
}





void extractRows(mat2 &A, mat2 &B, int start, int end) {
    int cantRows = end - start;
    for (int i = 0; i < cantRows; ++i) {
        for (int j = 0; j < B.cols(); ++j) {
            long double val = B.at(start + i, j);
            A.set(i, j, val);
        }
    }
}



void copyRowsTo(mat2 &A, mat2 &B, int start) {
    for (int i = 0; i < B.rows(); ++i) {
        for (int j = 0; j < B.cols(); ++j) {
            long double val = B.at(i, j);
            A.set(i + start, j, val);
        }
    }
}


void setPBorders(int myId, int cantProcs, mat2 &P0, mat2 &P1, mat2 &P2, int nX, int nY) {
    for (int r = 0; r < nX; ++r) {
        P0.set(r, 0, P0.at(r, 1));
        P1.set(r, 0, P1.at(r, 1));
        P2.set(r, 0, P2.at(r, 1));
    }

    for (int r = 0; r < nX; ++r) {
        P0.set(r, nY - 1, P0.at(r, nY - 2));
        P1.set(r, nY - 1, P1.at(r, nY - 2));
        P2.set(r, nY - 1, P2.at(r, nY - 2));
    }
    if (myId == cantProcs - 1) {
        for (int c = 0; c < nY; ++c) {
            P0.set(nX - 1, c, P0.at(nX - 2, c));
            P1.set(nX - 1, c, P1.at(nX - 2, c));
            P2.set(nX - 1, c, P2.at(nX - 2, c));
        }
    }

    if (myId == 0) {
        for (int c = 0; c < nY; ++c) {
            P0.set(0, c, P0.at(1, c));
            P1.set(0, c, P1.at(1, c));
            P2.set(0, c, P2.at(1, c));
        }
    }
}


void process(int myId, int cantProcs, MPI_Status stat, mat2 &U0, mat2 &U1, mat2 &U2, mat2 &V0, mat2 &V1, mat2 &V2, mat2 &P0, mat2 &P1, mat2 &P2, long double fanAngle, vector<pair<int, int>> rpt) {
    int startRow = rpt[myId].first;
    int endRow = rpt[myId].second;

    int nXlocal = U0.rows(); //replace nX and nY with local ones.

    //pre-allocating matrixes for sharing border data.
    mat2 U1ShareU(1, nY);
    mat2 U1ShareD(1, nY);
    mat2 V1ShareU(1, nY);
    mat2 V1ShareD(1, nY);
    mat2 P1ShareU(1, nY);
    mat2 P1ShareD(1, nY);

    unsigned int step = 0;

    for (long double t = 0.0; t < tMax; t = t + dt) {
        step++;
        //clearScreen();
        long double dFanAngle = fanTurns * 2 * pi / nT; //
        fanAngle += dFanAngle; //TODO: solo funciona en el 1er y tercer cuadrante.

        if (fanAngle > 2 * pi) fanAngle = 0;
        //cerr << "P" << myId << " " <<  100 * t / tMax << "%" << endl ;
        if (isnan(U1.at(3, 3))) {
            cerr << "ERROR: nan found" << endl;
            exit(EXIT_FAILURE);
        }

        for (int k = 0; k != 1; ++k) {


            setPBorders(myId, cantProcs, P0, P1, P2, nXlocal, nY);


            //process from start+1 to end-1
            for (int i = 1; i < nXlocal - 1; ++i) {
                for (int j = 1; j < nY - 1; ++j) {

                    long double U1x = (U1.at(i + 1, j) - U1.at(i - 1, j)) / (2.0 * dx);
                    long double U2x = (U2.at(i + 1, j) - U2.at(i - 1, j)) / (2.0 * dx);
                    long double U1y = (U1.at(i, j + 1) - U1.at(i, j - 1)) / (2.0 * dy);
                    long double U2y = (U2.at(i, j + 1) - U2.at(i, j - 1)) / (2.0 * dy);
                    long double U1xx = (U1.at(i + 1, j) - 2.0 * U1.at(i, j) + U1.at(i - 1, j)) / (dx * dx);
                    long double U2xx = (U2.at(i + 1, j) - 2.0 * U2.at(i, j) + U2.at(i - 1, j)) / (dx * dx);
                    long double U1yy = (U1.at(i, j + 1) - 2.0 * U1.at(i, j) + U1.at(i, j - 1)) / (dy * dy);
                    long double U2yy = (U2.at(i, j + 1) - 2.0 * U2.at(i, j) + U2.at(i, j - 1)) / (dy * dy);
                    long double P1x = (P1.at(i + 1, j) - P1.at(i - 1, j)) / (2.0 * dx);
                    long double P2x = (P2.at(i + 1, j) - P2.at(i - 1, j)) / (2.0 * dx);
                    long double P1y = (P1.at(i, j + 1) - P1.at(i, j - 1)) / (2.0 * dy);
                    long double P2y = (P2.at(i, j + 1) - P2.at(i, j - 1)) / (2.0 * dy);
                    long double V1x = (V1.at(i + 1, j) - V1.at(i - 1, j)) / (2.0 * dx);
                    long double V2x = (V2.at(i + 1, j) - V2.at(i - 1, j)) / (2.0 * dx);
                    long double V1y = (V1.at(i, j + 1) - V1.at(i, j - 1)) / (2.0 * dy);
                    long double V2y = (V2.at(i, j + 1) - V2.at(i, j - 1)) / (2.0 * dy);
                    long double V1xx = (V1.at(i + 1, j) - 2.0 * V1.at(i, j) + V1.at(i - 1, j)) / (dx * dx);
                    long double V2xx = (V2.at(i + 1, j) - 2.0 * V2.at(i, j) + V2.at(i - 1, j)) / (dx * dx);
                    long double V1yy = (V1.at(i, j + 1) - 2.0 * V1.at(i, j) + V1.at(i, j - 1)) / (dy * dy);
                    long double V2yy = (V2.at(i, j + 1) - 2.0 * V2.at(i, j) + V2.at(i, j - 1)) / (dy * dy);
                    long double P1xx = (P1.at(i + 1, j) - 2.0 * P1.at(i, j) + P1.at(i - 1, j)) / (dx * dx);
                    long double P1yy = (P1.at(i, j + 1) - 2.0 * P1.at(i, j) + P1.at(i, j - 1)) / (dy * dy);
                    long double P2yy = (P2.at(i, j + 1) - 2.0 * P2.at(i, j) + P2.at(i, j - 1)) / (dy * dy);



                    long double u2val = U1.at(i, j) - U1.at(i, j) * (dt / dx) * (U1.at(i, j) - U1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (U1.at(i, j) - U1.at(i, j - 1))
                                        - (dt / (rho * 2 * dx)) * (P1.at(i + 1, j) - P1.at(i - 1, j))
                                        + nu * ((dt / (dx * dx)) * (U1.at(i + 1, j) - 2 * U1.at(i, j) + U1.at(i - 1, j)) + (dt / (dy * dy)) * (U1.at(i, j + 1) - 2 * U1.at(i, j) + U1.at(i, j - 1)));
                    U2.set(i, j, u2val);

                    long double v2val = V1.at(i, j) - U1.at(i, j) * (dt / dx) * (V1.at(i, j) - V1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (V1.at(i, j) - V1.at(i, j - 1))
                                        - (dt / (rho * 2 * dy)) * (P1.at(i, j + 1) - P1.at(i, j - 1))
                                        + nu * ((dt / (dx * dx)) * (V1.at(i + 1, j) - 2 * V1.at(i, j) + V1.at(i - 1, j)) + (dt / (dy * dy)) * (V1.at(i, j + 1) - 2 * V1.at(i, j) + V1.at(i, j - 1)));
                    V2.set(i, j, v2val);

                    long double finalTerm = (1 / dt) * (U1x + V1y) - U1x * U1x - 2 * U1y * V1x - V1y * V1y;
                    long double p2val = ((P1.at(i + 1, j) + P1.at(i - 1, j)) * dy * dy + (P1.at(i, j + 1) + P1.at(i, j - 1)) * dx * dx) * (1.0 / (2 * (dx * dx + dy * dy)));
                    p2val -= (rho * dx * dx * dy * dy / (2 * dx * dx + 2 * dy * dy)) * finalTerm;
                    P2.set(i, j, p2val);

                    long double x = (i + startRow) * dx - xc;
                    long double y = j * dy - yc;
                    long double theta = atan(y / x);
                    long double beta = theta + (pi / 2.0);

                    long double r = sqrt(x * x + y * y);
                    long double tangSpeed = dFanAngle * r;

                    long double F = 2.0;
                    long double Fu = F * r * cos(beta);
                    long double Fv = F * r * sin(beta);

                    if (i < nXlocal / 2.0) {
                        Fu = -Fu;
                        Fv = -Fv;
                    }
                    if ( fabs(theta - fanAngle) < 0.1  && r > rMin && r < rMax) {
                        U2.add(i, j, Fu * dt);
                        V2.add(i, j, Fv * dt);
                    }


                }

            }


        }

        //Updating values
        U0.setAll(U1);
        U1.setAll(U2);
        V0.setAll(V1);
        V1.setAll(V2);
        P0.setAll(P1);
        P1.setAll(P2);


        //sharing borders. Borders are start+1 and end-1

        if (myId != 0) {
            extractRows(U1ShareU, U1, 1, 2);
            MPI_Send(U1ShareU.data, U1ShareU.cols(), MPI_LONG_DOUBLE, myId - 1, TAG, MPI_COMM_WORLD);
        }
        if (myId != cantProcs - 1) {
            extractRows(U1ShareD, U1, nXlocal - 2, nXlocal - 1);
            MPI_Send(U1ShareD.data, U1ShareD.cols(), MPI_LONG_DOUBLE, myId + 1, TAG, MPI_COMM_WORLD);
        }

        if (myId != 0) {

            extractRows(V1ShareU, V1, 1, 2);
            MPI_Send(V1ShareU.data, V1ShareU.cols(), MPI_LONG_DOUBLE, myId - 1, TAG, MPI_COMM_WORLD);
        }

        if (myId != cantProcs - 1) {

            extractRows(V1ShareD, V1, nXlocal - 2, nXlocal - 1);
            MPI_Send(V1ShareD.data, V1ShareD.cols(), MPI_LONG_DOUBLE, myId + 1, TAG, MPI_COMM_WORLD);
        }
        if (myId != 0) {

            extractRows(P1ShareU, P1, 1, 2);
            MPI_Send(P1ShareU.data, P1ShareU.cols(), MPI_LONG_DOUBLE, myId - 1, TAG, MPI_COMM_WORLD);
        }
        if (myId != cantProcs - 1) {

            extractRows(P1ShareD, P1, nXlocal - 2, nXlocal - 1);
            MPI_Send(P1ShareD.data, P1ShareD.cols(), MPI_LONG_DOUBLE, myId + 1, TAG, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //Get borders from neighbors, borders are start and end.
        long double buff[500];

        if (myId != 0) {
            MPI_Recv(U1ShareU.data, U1ShareU.cols(), MPI_LONG_DOUBLE, myId - 1, TAG, MPI_COMM_WORLD, &stat);
            U1.setRow(0, U1ShareU);
        }

        if (myId != cantProcs - 1) {
            MPI_Recv(U1ShareD.data, U1ShareD.cols(), MPI_LONG_DOUBLE, myId + 1, TAG, MPI_COMM_WORLD, &stat);
            U1.setRow(nXlocal - 1,  U1ShareD);
        }

        if (myId != 0) {

            MPI_Recv(V1ShareU.data, V1ShareU.cols(), MPI_LONG_DOUBLE, myId - 1, TAG, MPI_COMM_WORLD, &stat);
            V1.setRow(0, V1ShareU);
        }

        if (myId != cantProcs - 1) {

            MPI_Recv(V1ShareD.data, V1ShareD.cols(), MPI_LONG_DOUBLE, myId + 1, TAG, MPI_COMM_WORLD, &stat);
            V1.setRow(nXlocal - 1,  V1ShareD);
        }

        if (myId != 0) {

            MPI_Recv(P1ShareU.data, P1ShareU.cols(), MPI_LONG_DOUBLE, myId - 1, TAG, MPI_COMM_WORLD, &stat);
            P1.setRow(0, P1ShareU);
        }
        if (myId != cantProcs - 1) {

            MPI_Recv(P1ShareD.data, P1ShareD.cols(), MPI_LONG_DOUBLE, myId + 1, TAG, MPI_COMM_WORLD, &stat);
            P1.setRow(nXlocal - 1,  P1ShareD);
        }

        //si n%algo=0, enviar a rank 0 para imprimir.



        //si n%algo=0, enviar a rank 0 para imprimir.
        if (step % stepsUntilPrint == 0) {
            MPI_Barrier(MPI_COMM_WORLD);


            if (myId == 0) {
                mat2 U(nX, nY);
                copyRowsTo(U, U1, 0);

                mat2 V(nX, nY);
                copyRowsTo(V, U1, 0);

                for (int p = 1; p < cantProcs; ++p) {
                    int rowRange = rpt[p].second - rpt[p].first;
                    mat2 tmp(rowRange, nY);
                    MPI_Recv(tmp.data, rowRange * nY, MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD, &stat);
                    copyRowsTo(U, tmp, rpt[p].first);

                    MPI_Recv(tmp.data, rowRange * nY, MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD, &stat);
                    copyRowsTo(V, tmp, rpt[p].first);
                }
                U.print();
                V.print();

            } else {
                MPI_Send(U1.data, U1.rows()*U1.cols(), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD);
                MPI_Send(V1.data, V1.rows()*V1.cols(), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD);

            }

        }


    }
}



vector<pair<int, int>> calcRowsPerThread(int rows, int cantProcs) {
    int rowsPerThread = rows / (cantProcs - 1); //Watch out for the last one.
    //cout << "rows, procs, rowsPerThread = " << rows << ", " << cantProcs << ", " << rowsPerThread << endl;
    int remainingRows = rows - (cantProcs - 1) * rowsPerThread;
    vector<pair<int, int>> threadToRows;

    for (int i = 0; i < cantProcs; ++i) {
        int start = remainingRows + (i - 1) * rowsPerThread - 1;
        int end = remainingRows + i * rowsPerThread;
        if (i == 0) {
            start = 0;
            end = remainingRows;
        }

        if (i == cantProcs - 1) end--;
        cerr << start << ", " << end << endl;
        pair<int, int> rows = make_pair (start, end);
        threadToRows.push_back(rows);
    }
    return threadToRows;

}






int main(int argc, char *argv[]) {
    int cantProcs;
    int myId;
    MPI_Status stat;
    /* MPI programs start with MPI_Init; all 'N' processes exist thereafter */
    MPI_Init(&argc, &argv);
    /* find out how big the SPMD world is */
    MPI_Comm_size(MPI_COMM_WORLD, &cantProcs);
    /* and this processes' rank is */
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    /* At this point, all programs are running equivalently, the rank
       distinguishes the roles of the programs in the SPMD model, with
       rank 0 often used specially... */


    cerr << "nX, nY = " << nX << ", " << nY << endl;
    unsigned long int iter = 0;
    long double fanAngle = 0.0;
    vector<pair<int, int>> rpt = calcRowsPerThread(nX, cantProcs);


    if (myId == 0) {
        //RANK == 0.

        mat2 U0(nX, nY, 0.0);
        mat2 V0(nX, nY, 0.0);
        mat2 P0(nX, nY, 0.0);
        mat2 U1(nX, nY, 0.0);
        mat2 V1(nX, nY, 0.0);
        mat2 P1(nX, nY, 0.0);
        mat2 U2(nX, nY, 0.0);
        mat2 V2(nX, nY, 0.0);
        mat2 P2(nX, nY, 0.0);

        //Split and send:
        int cantRows;

        for (int p = 1; p < cantProcs; p++) {

            MPI_Send(&rpt[p].first, 1, MPI_INT, p, TAG, MPI_COMM_WORLD);
            MPI_Send(&rpt[p].second, 1, MPI_INT, p, TAG, MPI_COMM_WORLD);

            cantRows = rpt[p].second - rpt[p].first;

            mat2 U0frag(cantRows, nY);
            extractRows(U0frag, U0, rpt[p].first, rpt[p].second);
            MPI_Send(U0frag.data, U0frag.cols()*U0frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 U1frag(cantRows, nY);
            extractRows(U1frag, U1, rpt[p].first, rpt[p].second);
            MPI_Send(U1frag.data, U1frag.cols()*U1frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 U2frag(cantRows, nY);
            extractRows(U2frag, U2, rpt[p].first, rpt[p].second);
            MPI_Send(U2frag.data, U2frag.cols()*U2frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 V0frag(cantRows, nY);
            extractRows(V0frag, V0, rpt[p].first, rpt[p].second);
            MPI_Send(V0frag.data, V0frag.cols()*V0frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 V1frag(cantRows, nY);
            extractRows(V1frag, V1, rpt[p].first, rpt[p].second);
            MPI_Send(V1frag.data, V1frag.cols()*V1frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 V2frag(cantRows, nY);
            extractRows(V2frag, V2, rpt[p].first, rpt[p].second);
            MPI_Send(V2frag.data, V2frag.cols()*V2frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 P0frag(cantRows, nY);
            extractRows(P0frag, P0, rpt[p].first, rpt[p].second);
            MPI_Send(P0frag.data, P0frag.cols()*P0frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 P1frag(cantRows, nY);
            extractRows(P1frag, P1, rpt[p].first, rpt[p].second);
            MPI_Send(P1frag.data, P1frag.cols()*P1frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

            mat2 P2frag(cantRows, nY);
            extractRows(P2frag, P2, rpt[p].first, rpt[p].second);
            MPI_Send(P2frag.data, P2frag.cols()*P2frag.rows(), MPI_LONG_DOUBLE, p, TAG, MPI_COMM_WORLD);

        }

        cantRows = rpt[0].second - rpt[0].first;

        mat2 U0l(cantRows, nY, 0);
        mat2 U1l(cantRows, nY, 0);
        mat2 U2l(cantRows, nY, 0);
        mat2 V0l(cantRows, nY, 0);
        mat2 V1l(cantRows, nY, 0);
        mat2 V2l(cantRows, nY, 0);
        mat2 P0l(cantRows, nY, 0);
        mat2 P1l(cantRows, nY, 0);
        mat2 P2l(cantRows, nY, 0);

        extractRows(U0l, U0, rpt[0].first, rpt[0].second);
        extractRows(U1l, U1, rpt[0].first, rpt[0].second);
        extractRows(U2l, U2, rpt[0].first, rpt[0].second);
        extractRows(V0l, V0, rpt[0].first, rpt[0].second);
        extractRows(V1l, V1, rpt[0].first, rpt[0].second);
        extractRows(V2l, V2, rpt[0].first, rpt[0].second);
        extractRows(P0l, P0, rpt[0].first, rpt[0].second);
        extractRows(P1l, P1, rpt[0].first, rpt[0].second);
        extractRows(P2l, P2, rpt[0].first, rpt[0].second);

        process(myId, cantProcs, stat,  U0l, U1l, U2l, V0l, V1l, V2l, P0l, P1l, P2l, fanAngle, rpt);



    } else {
        //RANK != 0.

        //recieve:
        int startRow;
        int endRow;
        MPI_Recv(&startRow, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat);
        MPI_Recv(&endRow, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat);


        int cantRows = endRow - startRow;

        cerr << "p " << myId << ", cantRows " << cantRows << endl;

        mat2 U0l(cantRows, nY, 0);
        MPI_Recv(U0l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 U1l(cantRows, nY, 0);
        MPI_Recv(U1l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 U2l(cantRows, nY, 0);
        MPI_Recv(U2l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 V0l(cantRows, nY, 0);
        MPI_Recv(V0l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 V1l(cantRows, nY, 0);
        MPI_Recv(V1l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 V2l(cantRows, nY, 0);
        MPI_Recv(V2l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 P0l(cantRows, nY, 0);
        MPI_Recv(P0l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 P1l(cantRows, nY, 0);
        MPI_Recv(P1l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        mat2 P2l(cantRows, nY, 0);
        MPI_Recv(P2l.data, nY * (endRow - startRow), MPI_LONG_DOUBLE, 0, TAG, MPI_COMM_WORLD, &stat);

        process(myId, cantProcs, stat,  U0l, U1l, U2l, V0l, V1l, V2l, P0l, P1l, P2l, fanAngle, rpt);
    }

    MPI_Finalize();
    return 0;
}
