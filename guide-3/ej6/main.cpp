#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include <math.h>


using namespace std;

using mat3 = vector<vector<vector<long double > > >;
using mat2 = vector<vector<long double > >;
using mat1 = vector<long double>;


int main() {
    //Donde x, y ∈ [0; 2] y t > 0 son las variables espacial y temporal respectivamente. u(x, y, t) y v(x, y, t) son los
// componentes x-y de la velocidad del fluido dependiente del espacio y el tiempo. Se desea resolver el Lid −
// driven cavity problem en el cual a un recipiente (2D) que contiene un determinado fluido se le desliza su tapa
// hacia uno de los lados, provocando cambios en el movimiento y la presi ́on del fluido. Para resolver el problema
// se cuenta con los siguiente datos: ν = 0.1. ρ=1. Condiciones iniciales u, v, p = 0 en todos lados. Condiciones
// ∂p
// ∂p
// de contorno u = 1 en y = 2 (la “tapa”). u, v = 0 en los dem ́as bordes. ∂y
// = 0 en y = 0. p = 0 en y = 2. ∂x
// =0
// en x = 0, 2
    long double xMax = 2.0;
    long double yMax = 2.0;
    long double tMax = 1.0;
    long double nu = 0.1; //viscosidad
    long double rho = 1;   //densidad
    long double dx = 0.2;
    long double dy = 0.2;
    long double dt = 0.5;
    int nX = round(xMax / dx) + 1;
    int nY = round(yMax / dy) + 1;
    int nT = round(tMax / dt) + 1;
    long double al = 0.5;
    mat2 U0(nX, vector<long double>(nY));
    mat2 V0(nX, vector<long double>(nY));
    mat2 P0(nX, vector<long double>(nY));
    mat2 U1(nX, vector<long double>(nY));
    mat2 V1(nX, vector<long double>(nY));
    mat2 P1(nX, vector<long double>(nY));
    mat2 U2(nX, vector<long double>(nY));
    mat2 V2(nX, vector<long double>(nY));
    mat2 P2(nX, vector<long double>(nY));
    for (int i = 0; i < nX; ++i) {
        U0[i][nY-1] = 1.0;
        U1[i][nY-1] = 1.0;
        U2[i][nY-1] = 1.0;
    }
    for (int i = 1; i < nX - 1; ++i) {
        for (int j = 1; i < nY - 1; ++i) {
            for (int k = 0;    k < 10; ++k) {
                long double U1x = (U1[i + 1][j] - U1[i - 1][j]) / (2 * dx);
                long double U2x = (U2[i + 1][j] - U2[i - 1][j]) / (2 * dx);
                long double U1y = (U1[i][j + 1] - U1[i][j - 1]) / (2 * dy);
                long double U2y = (U2[i][j + 1] - U2[i][j - 1]) / (2 * dy);
                long double U1xx = (U1[i + 1][j] - 2 * U1[i][j] + U1[i - 1][j]) / (dx * dx);
                long double U2xx = (U2[i + 1][j] - 2 * U2[i][j] + U2[i - 1][j]) / (dx * dx);
                long double U1yy = (U1[i][j + 1] - 2 * U1[i][j] + U1[i][j - 1]) / (dy * dy);
                long double U2yy = (U2[i][j + 1] - 2 * U2[i][j] + U2[i][j - 1]) / (dy * dy);
                long double P1x = (P1[i + 1][j] - P1[i - 1][j]) / (2 * dx);
                long double P2x = (P2[i + 1][j] - P2[i - 1][j]) / (2 * dx);
                long double P1y = (P1[i][j + 1] - P1[i][j - 1]) / (2 * dy);
                long double P2y = (P2[i][j + 1] - P2[i][j - 1]) / (2 * dy);
                long double V1x = (V1[i + 1][j] - V1[i - 1][j]) / (2 * dx);
                long double V2x = (V2[i + 1][j] - V2[i - 1][j]) / (2 * dx);
                long double V1y = (V1[i][j + 1] - V1[i][j - 1]) / (2 * dy);
                long double V2y = (V2[i][j + 1] - V2[i][j - 1]) / (2 * dy);
                long double V1xx = (V1[i + 1][j] -    2 * V1[i][j] + V1[i - 1][j]) / (dx * dx);
                long double V2xx = (V2[i + 1][j] - 2 * V2[i][j] + V2[i - 1][j]) / (dx * dx);
                long double V1yy = (V1[i][j + 1] - 2 * V1[i][j]    + V1[i][j - 1]) / (dy * dy);
                long double V2yy = (V2[i][j + 1] - 2 * V2[i][j] + V2[i][j - 1]) / (dy * dy);
                long double P1xx = (P1[i + 1][j] - 2 * P1[i][j] + P1[i - 1][j]) / (dx * dx);
                long double P1yy = (P1[i][j + 1] - 2 * P1[i][j] + P1[i][j - 1]) / (dy * dy);
                long double P2yy = (P2[i][j + 1] - 2 * P2[i][j] + P2[i][j - 1]) / (dy * dy);
                //TODO: chequear despejes de ecuación 3
                U2[i][j] = U0[i][j] + 2 * dt * (-U1[i][j] * (al * U1x  + (1 - al) * U2x ) - V1[i][j] * (al * U1y + (1 - al) * U2y)
                                                - (1 / rho) * (al * P1x + (1 - al) * P2x) + nu * (al * U1xx + (1 - al) * U2xx + al * U1yy + (1 - al) * U2yy));
                V2[i][j] = V0[i][j] + 2 * dt * (-U1[i][j] * (al * V1x  + (1 - al) * V2x ) - V1[i][j] * (al * V1y + (1 - al) * V2y)
                                                - (1 / rho) * (al * P1y + (1 - al) * P2y) + nu * (al * V1xx + (1 - al) * V2xx + al * V1yy + (1 - al) * V2yy));
                P2[i][j] = (P2[i - 1][j] / 2)  + (P2[i + 1][j] / 2) + (dx * dx) / (2 * al - 2) * (-al * P1xx - al * P1yy - (1 - al) * P2yy
                           - rho * ( pow(2, (al * U1x + (1 - al) * U2x)) + 2 * (al * U1y + (1 - al) * U2y) * (al * V1x + (1 - al) * V2x) + pow(2, (al * V1y + (1 - al) * V2y))));

            }
        }
    }
}



