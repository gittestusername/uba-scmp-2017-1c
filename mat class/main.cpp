#include <iostream>
#include <vector>
#include "mat.h"

using namespace std;

int main() {

    mat<double> A(5, 5, 1);
    A.id();

    mat<double> B(5, 1, 0);
    B.id();
    B.set(0, 1);
    B.set(1, 2);
    B.set(2, 3);




    mat<double> x1 = A.jacobi(B);
    cout << x1;
    mat<double> x2 = A.gaussElimination(B);
    cout << x2;
    return 0;

}