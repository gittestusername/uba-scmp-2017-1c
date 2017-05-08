#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include "../../mat/mat.h"

using namespace std;

int main() {
    /*
        long double m = 1.0;             //kg.
        long double l = 1.0;             //meters.
        long double mu = m / l;
        long double tension = 40.0;      //kgf. Tension.
        long double k = tension / mu;
        long double pullPoint = 0.5;    //meters from left.
        long double pullLength = 0.07;   //pull meters.
        long double xMax = 1.0;
        long double tMax = 1.0;          //tMin and xMin are zero.
        long double dx = 0.005;
        long double dt = 0.3;
    */


    long double m = 1.0 / 1000.0;         //kg.
    long double xMax = 0.6;            //meters
    long double tMax = 0.6;            //tMin and xMin are zero.
    long double mu = m / xMax;
    long double tension = 40.0 * 9.80665;    //Newton. Tension.
    long double k = tension / mu;
    long double pullPoint = 0.5;    //meters from left.
    long double pullLength = 0.07;   //pull meters.

    long double dx = 0.2;
    long double dt = 0.2;

    int nX = round(xMax / dx) + 1; //if dx = 0.5, we need 0,0.5, and 1.0 values.
    int nT = round(tMax / dt) + 1;
    int n = nX * nT;





    mat<long double> A(n, n);
    mat<long double> B(n, 1, 0);
    int schemeSize = 4 * nX + 1;
    int lenTails = (schemeSize - 1) / 2;


    mat<long double> scheme(1, schemeSize, 0);



    /*
        scheme.set(lenTails, 2.0 - 2.0 * (dt * dt) / (k * dx * dx));
        scheme.set(lenTails - 1.0, (dt * dt) / (k * dx * dx));
        scheme.set(lenTails + 1.0, (dt * dt) / (k * dx * dx));
        scheme.set(lenTails - nX, -1.0);
    */
//From class.



    scheme.set(lenTails, 2.0 * (1 - k * (dt * dt) / (dx * dx)));
    scheme.set(lenTails - 1, k * (dt * dt) / (dx * dx));
    scheme.set(lenTails + 1, k * (dt * dt) / (dx * dx));
    scheme.set(lenTails - nX, -1.0);


    /*From class,  simplified (asumes dt = dx = 1.)
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
        long double slope = pullLength / pullPoint;
        long double u;
        if (x < pullPoint) {
            u = slope * x;
        } else {
            u = pullLength - slope * (x - pullPoint);
        }
        if(abs(u) < 0.00001) u = 0.0;
        B.set(i, u);
        cout << " X, u = " << x  << ", " << u << endl;
        x += dx;
    }

cout << B << endl;
//We need now to set the A and B so that x0 and xf ar always 0.

    for (int i = 1; i < nT; ++i) {
        int row = i * nX;
        A.setRow(row, 0);
        A.set(row, row, 1);

        A.setRow(row - 1, 0);
        A.set(row - 1, row - 1, 1);


        B.set(row, 0);
        //B.set(row-1, 0);

    }
    A.setRow(n - 1, 0);
    A.set(n - 1, n - 1, 1);
    B.set(n - 1, 0);


    cout << A << endl;

    mat<long double> X = A.gaussElimination(B);
    int count = 0;
    cout << " X ROWS = n?, " << X.rows() <<  " = " << n << endl;
    cout << " nX, nT =  " << nX <<  ", " << nT << endl;
    cout << B << endl;
    for (int i = 0; i < X.rows(); ++i) {
        count++;
        if(abs(X.at(i)) < 0.00001) X.set(i,0);
        cout << X.at(i) << " ";

        if (count == nX && i < X.rows() - 1) {
            count = 0;
            cout << endl;
        }
    }



}