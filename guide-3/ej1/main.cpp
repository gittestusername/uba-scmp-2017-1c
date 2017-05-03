#include <iostream>
#include <vector>
#include "../../mat/mat.h"

using namespace std;

int main() {

//Differential equation is d2u/dt2 = (T/'mu') * d2u/dx2
// where x in [0;1] and t > 0 are the spatial and temporal variables.
// u(x,t)


    double delta = 1;
    double max = 10;
    int n = max / delta;

    mat<double> A(n, n, 0);
    //its n, not n*n, because its unidimensional.
    mat<double> scheme(1, 3, 0);
    scheme.set(0, 1 , 1);
    scheme.set(0, 2 , 1.0 / 4.0);
    scheme.set(0, 0 , 1.0 / 4.0);

    A.fillScheme(scheme);
    cout << A << endl;

    return 0;

}