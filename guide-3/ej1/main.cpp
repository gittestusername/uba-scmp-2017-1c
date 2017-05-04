#include <iostream>
#include <vector>
#include "../../mat/mat.h"

using namespace std;

int main() {

//Differential equation is d2u/dt2 = (T/'mu') * d2u/dx2
// where x in [0;1] and t > 0 are the spatial and temporal variables.
// u(x,t)


    double delta = 1;
    double max = 3;
    int n = max / delta;

    mat<double> A(n, 2*n, 3);
    mat<double> B(2*n, n, 3);

    int tmp = 1;
for (int i = 0; i < A.rows(); ++i)
{
    for (int j = 0; j < A.cols(); ++j)
    {
        A.set(i,j, tmp);
        B.set(j,i, (tmp)%3 + tmp);
        tmp++;
    }
}

    mat<double> C1 = (A.sparseProd(B));
    mat<double> C2 = A*B;
    mat<double> dif = C1-C2;


    cout << "A" << endl << endl << A << endl;

    cout << "B" << endl << endl << B << endl;
    cout << "C1" << endl << endl << C1 << endl;

    cout << "C2" << endl << endl << C2 << endl;
    cout << "dif" << endl << endl << dif << endl;

    return 0;

}