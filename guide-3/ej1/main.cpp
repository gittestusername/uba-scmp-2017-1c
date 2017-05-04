#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock
#include "../../mat/mat.h"

using namespace std;

int main() {

//Differential equation is d2u/dt2 = (T/'mu') * d2u/dx2
// where x in [0;1] and t > 0 are the spatial and temporal variables.
// u(x,t)

    int n = 3000;
    mat<double> scheme(1,201,1);


    mat<double> A(n, n, 0);
    mat<double> B(n, 1, 1);

    A.fillScheme(scheme);


    auto t1 = std::chrono::high_resolution_clock::now();
    mat<double> C1  = A*B;
    auto t2 = std::chrono::high_resolution_clock::now();
 
    // integral duration: requires duration_cast
    auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
 
    // fractional duration: no duration_cast needed
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
 
    std::cout << "regular took " << int_ms.count() << " whole milliseconds\n";








    t1 = std::chrono::high_resolution_clock::now();
    mat<double> C2  = A.sparseProd(B);
    t2 = std::chrono::high_resolution_clock::now();
 
    // integral duration: requires duration_cast
    int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
 
    // fractional duration: no duration_cast needed
    fp_ms = t2 - t1;
 
    std::cout << "sparse took " << int_ms.count() << " whole milliseconds\n";



    return 0;

}