#include <iostream>
#include <vector>
#include "mat.h"

using namespace std;

int main() {

    double delta = 1;
    double max = 10;
    int n = max/delta;

    mat<double> A(n, n, 0);
    //its n, not n*n, because its unidimensional.
    mat<double> scheme(1, n, 0);
    scheme.set(0,n/2, 1);
    scheme.set(0,n/2 + 1 , 1/4);
    scheme.set(0,n/2 - 1 , 1/4);
    cout << scheme.at(0,5) << endl;


    cout << scheme << endl;
    return 0;

}