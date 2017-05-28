#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include "../../../mat/mat.h"

using namespace std;

int main() {
    long double m = 1.0 / 1000.0;             //kg.
    long double xMax = 1.0;                     //meters
    long double tMax = 0.2;                     //tMin and xMin are zero.
    long double mu = m / xMax;
    long double tension = 40.0*9.80665;     //Newton. Tension.
    long double pullX = 0.5;                  //meters from left.
    long double pull = 0.07;                  //pull meters.
    long double dx = 0.05;
    long double dt = 0.000001;
    //For reference, dx = 0.005 y dt = 0.0000002 doesnt finish, too long
    int nX = round(xMax / dx) + 1;
    int nT = round(tMax / dt) + 1;
    //if dx = 0.5, we need the values of x = 0,
    // x = 0.5, and x = 1.0. so there is one more.
    //the same happens for time.
    long double s = (tension / mu)*(dt*dt) / ( dx*dx);
    mat<long double> U0(nX, 1, 0);
    mat<long double> U1(nX, 1, 0);
    long double x = 0;
    for (int i = 0; i < nX; ++i) {
        //Fill U0 and U1
        long double slope = pull / pullX;
        long double u;
        if (x < pullX) {
            u = slope*x;
        } else {
            u = pull - slope*(x - pullX);
        }
        U0.set(i, u);
        U1.set(i, u);
        x += dx;
    }
    //Defining the scheme. 
    long double alpha = 0.5;
    long double k = (mu*dx*dx) / (tension*dt*dt);
    //fill A
    mat<long double> schemeA(1, 3, 0);
    schemeA.set(0, alpha - 1);
    schemeA.set(1, k + 2*(1-alpha) );
    schemeA.set(2, alpha - 1);
    mat<long double> A(nX, nX, 0);
    A.fillScheme(schemeA);
    //fill B
    mat<long double> schemeB(1, 3, 0);
    schemeB.set(0, alpha);
    schemeB.set(1, 2*k - 2*alpha);
    schemeB.set(2, alpha);
    mat<long double> B(nX, nX, 0);
    B.fillScheme(schemeB);
    //fill C
    mat<long double> schemeC(1, 1, 0);
    schemeC.set(0, -k);
    mat<long double> C(nX, nX, 0);
    C.fillScheme(schemeC);
    //Now i have to specify the border conditions
    A.setRow(0, 0);
    A.setRow(nX - 1, 0);
    A.set(0, 0, 1);
    A.set(nX - 1, nX - 1, 1);
    B.setRow(0, 0);
    B.setRow(nX - 1, 0);
    B.set(0, 0, 1);
    B.set(nX - 1, nX - 1, 1);
    C.setRow(0, 0);
    C.setRow(nX - 1, 0);

    for (int t = 0; t < nT; ++t) {
        mat<long double> BU = B*U1;
        mat<long double> CU = C*U0;
        mat<long double> K = BU + CU;
        mat<long double> R = A.gaussElimination(K);
        U0 = U1;
        U1 = R;
        //print result:
        cout << endl;
        for (int i = 0; i < U0.rows(); ++i) {
            cout << " " << U0.at(i);
        }
    }
}