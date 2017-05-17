#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include "../../mat/mat.h"

using namespace std;

int main() {

    long double m = 1.0 / 1000.0;             //kg.
    long double xMax = 1.0;                     //meters
    long double tMax = 0.01;                     //tMin and xMin are zero.
    long double mu = m / xMax;
    long double tension = 40.0 * 9.80665;     //Newton. Tension.
    long double pullX = 0.5;                  //meters from left.
    long double pull = 0.07;                  //pull meters.

    long double dx = 0.01;
    long double dt = 0.000001;

    int nX = round(xMax / dx) + 1;
    int nT = round(tMax / dt) + 1;
    //if dx = 0.5, we need the values of x = 0,
    // x = 0.5, and x = 1.0. so there is one more.
    //the same happens for time.
    long double s = (tension / mu)*(dt * dt) / ( dx * dx);

/*
    cout << "s = " << s << endl;
    cout << "tension = " << tension << endl;
    cout << "mu = " << mu << endl;
    cout << "m = " << m << endl;
    return 0;
*/

    mat<long double> U0(nX, 1, 0);
    mat<long double> U1(nX, 1, 0);


    long double x = 0;
    for (int i = 0; i < nX; ++i) {
 

        //Fill U0.
        long double slope = pull / pullX;
        long double u;
        if (x < pullX) {
            u = slope * x;
        } else {
            u = pull - slope * (x - pullX);
        }
        U0.set(i, u);
        U1.set(i, u);

        x += dx;
    }

    mat<long double> U2(nX, 1, 0);

    //cout << U0 << endl;

    for (int t = 0; t < nT; ++t) {
        if(true || t%20 == 3){
           for (int j = 0; j < U0.rows(); ++j)
           {
               cout << U0.at(j) << " ";
           }
           cout << endl;
        }
       for (int i = 1; i < nX-1; ++i)
       {
           //U2.set(i,s*(U1.at(i+1) -2*U1.at(i) + U1.at(i-1)) + 2*U1.at(i) + U0.at(i));
           U2.set(i, s*(U1.at(i+1) + U1.at(i-1)) + 2.0*(1.0-s)*(U1.at(i)) - U0.at(i));
       }
       U0 = U1;
       U1 = U2;
              //Calculate the period of the wave.
       long double max = 0.0;

      

/*
       for (int i = 0; i < U0.rows(); ++i)
       {
        if(U0.at(i) > max) max = U0.at(i);
       }


       if (t > 50 && max >= pull - 0.005)
       {
            cout << t*dt << endl;
            cout << " at iteration " <<  t << endl;
           return 0;
       }
*/



    }

    //cout << s << endl;
}