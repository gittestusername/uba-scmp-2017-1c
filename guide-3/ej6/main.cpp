#include <iostream>
#include <vector>
#include <chrono> // for high_resolution_clock.
#include <math.h>
#include <sstream>
#include "mat2.h"
#include "parameters.h"

using namespace std;


void clearScreen() {
    cerr << string( 100, '\n' );
}

long double diff(mat2 &A, mat2 &B) {
    long double sum = 0;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            sum += pow(B.at(i,j) - A.at(i,j), 2);
        }
    }
    return sum;
}

void printMat(mat2 m) {
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            cout << m.at(i,j) << " " ;
        }
        cout << endl;
    }
    cout << endl;
}

void setPBorders(mat2 &P0, mat2 &P1, mat2 &P2, int nX, int nY) {
    for (int c = 0; c < nX; ++c) {
        P0.set(c,0, P0.at(c,1));
        P1.set(c,0, P1.at(c,1));
        P2.set(c,0, P2.at(c,1));
    }
    for (int c = 0; c < nY; ++c) {
        P0.set(nX - 1,c, P0.at(nX - 2,c));
        P1.set(nX - 1,c, P1.at(nX - 2,c));
        P2.set(nX - 1,c, P2.at(nX - 2,c));
        P0.set(0,c, P0.at(1,c));
        P1.set(0,c, P1.at(1,c));
        P2.set(0,c, P2.at(1,c));
    }
}
// TODO: poner todo condiciones de Neumann, excepto en un solo punto del borde donde la presion nula.
//Esto es comp ara tener un cero de referencia.

int main() {

    //0 es tiempo n-1, 1 es tiempo n, 2 es tiempo n+1
    //3 es simplemente las matrices auxiliares para no perder datos de el tiempo 2
    mat2 U0(nX, nY, 0);
    mat2 V0(nX, nY, 0);
    mat2 P0(nX, nY, 0);
    mat2 U1(nX, nY, 0);
    mat2 V1(nX, nY, 0);
    mat2 P1(nX, nY, 0);
    mat2 U2(nX, nY, 0);
    mat2 V2(nX, nY, 0);
    mat2 P2(nX, nY, 0);
    mat2 U2_old(nX, nY, 0);
    mat2 V2_old(nX, nY, 0);
    mat2 P2_old(nX, nY, 0);
    mat2 tmp(nX, nY, 0);

    //condicion borde u=1 en y=2
    for (int i = 0; i < nX; ++i) {
        U0.set(i, nY - 1, 1.0);
        U1.set(i, nY - 1, 1.0);
        U2.set(i, nY - 1, 1.0);
    }
    for (long double t = 0.0; t < tMax; t = t + dt) {
        cerr << 100 * t / tMax << "%" << endl ;
        if (isnan(U1.at(3, 3))) {
            cerr << "ERROR: nan found" << endl;
            exit(EXIT_FAILURE);
        }
        U0.print();
        V0.print();
        //P0.print();
        tmp.setAll(U1);
        for (int k = 0; true; ++k) {
            tmp.setAll(U2);
            setPBorders(P0, P1, P2, nX, nY);

            //U2_old.setAll(U2);
            //V2_old.setAll(V2);
            //P2_old.setAll(P2);
            for (int i = 1; i < nX - 1; ++i) {
                for (int j = 1; j < nY - 1; ++j) { //las condiciones borde en Y=2 se respetan aca
                    //0 = n-1, 1 = n, 2 = n+1
                    //x = dx, y = dy, xx = d2x, yy = d2y.
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

                    //TODO: U0, V0, P0, are not used if we use forward scheme in time.
                    /*long double oU1x = U1x;
                    long double oU2x = U2x;
                    long double oV1x = V1x;
                    long double oV2x = V2x;
                    long double oU1y = U1y;
                    long double oU2y = U2y;
                    long double oV1y = V1y;
                    long double oV2y = V2y;*/

                    if (upwind) {
                        if (U1.at(i,j) > 0) {
                            U1x = (U1.at(i + 1, j) - U1.at(i, j)) / dt;
                            V1x = (V1.at(i + 1, j) - V1.at(i, j)) / dt;
                            U2x = (U2.at(i + 1, j) - U2.at(i, j)) / dt;
                            V2x = (V2.at(i + 1, j) - V2.at(i, j)) / dt;
                        } else {
                            U1x = (U1.at(i, j) - U1.at(i - 1, j)) / dt;
                            V1x = (V1.at(i, j) - V1.at(i - 1, j)) / dt;
                            U2x = (U2.at(i, j) - U2.at(i - 1, j)) / dt;
                            V2x = (V2.at(i, j) - V2.at(i - 1, j)) / dt;
                        }
                        if (V1.at(i,j) > 0) {
                            U1y = (V1.at(i, j + 1) - V1.at(i, j)) / dt;
                            V1y = (V1.at(i, j + 1) - V1.at(i, j)) / dt;
                            U2y = (V2.at(i, j + 1) - V2.at(i, j)) / dt;
                            V2y = (V2.at(i, j + 1) - V2.at(i, j)) / dt;
                        } else {
                            U1y = (V1.at(i, j) - V1.at(i, j - 1)) / dt;
                            V1y = (V1.at(i, j) - V1.at(i, j - 1)) / dt;
                            U2y = (V2.at(i, j) - V2.at(i, j - 1)) / dt;
                            V2y = (V2.at(i, j) - V2.at(i, j - 1)) / dt;
                        }
                    }


                    //Lorena Barba's discretization.
                    long double u2val = U1.at(i, j) - U1.at(i, j) * (dt / dx) * (U1.at(i, j) - U1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (U1.at(i, j) - U1.at(i, j - 1))
                                        - (dt / (rho * 2 * dx)) * (P1.at(i + 1, j) - P1.at(i - 1, j))
                                        + nu * ((dt / (dx * dx)) * (U1.at(i + 1, j) - 2 * U1.at(i, j) + U1.at(i - 1, j)) + (dt / (dy * dy)) * (U1.at(i, j + 1) - 2 * U1.at(i, j) + U1.at(i, j - 1)));
                    U2.set(i, j, u2val);

                    long double v2val = V1.at(i, j) - U1.at(i, j) * (dt / dx) * (V1.at(i, j) - V1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (V1.at(i, j) - V1.at(i, j - 1))
                                        - (dt / (rho * 2 * dy)) * (P1.at(i, j + 1) - P1.at(i, j - 1))
                                        + nu * ((dt / (dx * dx)) * (V1.at(i + 1, j) - 2 * V1.at(i, j) + V1.at(i - 1, j)) + (dt / (dy * dy)) * (V1.at(i, j + 1) - 2 * V1.at(i, j) + V1.at(i, j - 1)));
                    V2.set(i, j, v2val);

                    long double lastTerm = (1 / dt) * (U1x + V1y) - U1x * U1x - 2 * U1y * V1x - V1y * V1y;

                    long double p2val = ((P1.at(i + 1, j) + P1.at(i - 1, j)) * dy * dy + (P1.at(i, j + 1) + P1.at(i, j - 1)) * dx * dx) * (1.0 / (2 * (dx * dx + dy * dy)));
                    p2val -= (rho * dx * dx * dy * dy / (2 * dx * dx + 2 * dy * dy)) * lastTerm;
                    P2.set(i, j, p2val);

                }
            }

            long double matDiff = diff(tmp, U2);
            if (matDiff < fixedPointError && k > 5) {
                if (debug) cerr << "out at k = " << k << endl;
                break;
            } else if (k > 4000) {
                cerr << "ERROR: unstable." << endl;
                exit(EXIT_FAILURE);
            }
        }

        U0.setAll(U1);
        U1.setAll(U2);
        V0.setAll(V1);
        V1.setAll(V2);
        P0.setAll(P1);
        P1.setAll(P2);
    }
}
