#include <iostream>
#include <vector>
#include "mat.h"

using namespace std;

int main() {

    mat<double> A(2, 2, 0);
    mat<double> B(2, 2, 0);
    A.set(0, 0, 1);
    A.set(0, 1, 2);
    A.set(1, 0, 3);
    A.set(1, 1, 4);

    B.set(0, 0, 5);
    B.set(0, 1, 6);
    B.set(1, 0, 7);
    B.set(1, 1, 8);

    A = A * B;
    cout << A;

    return 0;

}