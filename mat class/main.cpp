#include <iostream>
#include <vector>
#include "mat.h"

using namespace std;

int main() {

    mat<double> A(20, 20, 0);
    mat<double> B(20, 20, 0);
    A.set(0,0,1.234);
    A.set(0,1,234);
    A.set(1,0,3);
    A.set(1,1,4);

	B.set(0,0,5);
    B.set(0,1,6);
    B.set(1,0,7);
    B.set(1,1,8);

    A = A.prod(B);
    A.print();
    return 0;

}