#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include "../../mat/mat.h"

using namespace std;

int main() {



    long double m = 1.0 / 1000.0;         //kg.
    long double xMax = 1;            //meters
    long double tMax = 3;            //tMin and xMin are zero.
    long double mu = m / xMax;
    long double tension = 40.0 * 9.80665;    //Newton. Tension.
    long double k = tension / mu;
    long double pullX = 0.5;    //meters from left.
    long double pull = 0.07;   //pull meters.

    long double dx = 0.1;
    long double dt = 0.2;

    int nX = round(xMax / dx) + 1;
    int nT = round(tMax / dt) + 1;
    //if dx = 0.5, we need the values of x = 0,
    // x = 0.5, and x = 1.0. so there is one more.
    //the same happens for time.
    int n = nX * nT;
    long double s = (dt * dt) / (k * dx * dx);




    mat<long double> A(n, n);
    mat<long double> B(n, 1, 0);
    int schemeSize = 4 * nX + 1;
    int lenTails = (schemeSize - 1) / 2;


    mat<long double> scheme(1, schemeSize, 0);

//My scheme
    scheme.set(lenTails, 2.0 - 2.0 * s);
    scheme.set(lenTails - 1, s);
    scheme.set(lenTails + 1, s);
    scheme.set(lenTails - nX, -1.0);

//From class.
    /*
       scheme.set(lenTails, 2.0 * (1 - k * (dt * dt) / (dx * dx)));
       scheme.set(lenTails - 1, k * (dt * dt) / (dx * dx));
       scheme.set(lenTails + 1, k * (dt * dt) / (dx * dx));
       scheme.set(lenTails - nX, -1.0);
    */


//From class,  simplified (asumes dt = dx = 1.)
    /*
        scheme.set(lenTails - 1, 1);
        scheme.set(lenTails + 1, 1);
        scheme.set(lenTails - nX, -1);
    */


    A.fillScheme(scheme);


    long double x = 0;
    for (int i = 0; i < nX; ++i) {
        A.setRow(i, 0);
        A.set(i, i, 1);

        //Fill B.
        long double slope = pull / pullX;
        long double u;
        if (x < pullX) {
            u = slope * x;
        } else {
            u = pull - slope * (x - pullX);
        }
        B.set(i, u);

        //We set a second point, because the scheme uses the last two iterations.
        int nxt = i + nX;
        A.setRow(nxt, 0);
        A.set(nxt, nxt, 1);
        B.set(nxt, u);

        x += dx;
    }

    //We need now to set the A and B so that x0 and xf ar always 0.
    for (int i = 1; i < nT; ++i) {
        int row = i * nX;
        A.setRow(row, 0);
        A.set(row, row, 1);

        A.setRow(row - 1, 0);
        A.set(row - 1, row - 1, 1);

        B.set(row, 0.00);
        B.set(row - 1, 0.00);

    }
    A.setRow(n - 1, 0.0);
    A.set(n - 1, n - 1, 1.0);
    B.set(n - 1, 0.0);

    //cout << B << endl;
    // cout << A << endl;

    mat<long double> X = A.gaussElimination(B);
    int count = 0;


    for (int i = 0; i < X.rows(); ++i) {
        count++;
        cout << X.at(i) << " ";

        if (count == nX && i < X.rows() - 1) {
            count = 0;
            cout << endl;
        }
    }



}