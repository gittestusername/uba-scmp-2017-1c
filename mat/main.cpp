#include "mat.h"

using namespace std;

int main() {

    mat<double> A(4, 4);
    A.id();
    A.set(0, 0, 10);
    A.set(0, 1, 5);
    A.set(0, 2, 7);
    A.set(0, 3, 4);
    A.set(1, 0, 5);
    A.set(1, 1, 15);
    A.set(1, 2, 7);
    A.set(1, 3, 3);
    A.set(2, 0, 9);
    A.set(2, 1, 3);
    A.set(2, 2, 11);
    A.set(2, 3, 4);
    A.set(3, 0, 1);
    A.set(3, 1, 2);
    A.set(3, 2, 5);
    A.set(3, 3, 14);


    //cout << A << endl;

    mat<double> Ai = A.inverse();
    mat<double> id = A * Ai;

    cout << id << endl;



    return 0;
}