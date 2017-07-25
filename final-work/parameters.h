using namespace std;

long double xMax = 15.0;
long double yMax = 15.0;
long double tMax = 15.0;
long double nu = 0.01; //viscosidad. Era 0.1.
long double rho = 1.0;  //densidad
long double dx = (1.0 / 20.0);
long double dy = (1.0 / 20.0);
long double dt = 0.001;
int nX = round(xMax / dx) + 1;
int nY = round(yMax / dy) + 1;
int nT = round(tMax / dt) + 1;
long double al = 0.5;
bool upwind = false;
long double fixedPointError = 0.000001;
long double minFixedPointIters = 10;
long double xc = xMax / 2;
long double yc = yMax / 2;
long double rMax = 0.5 * min(xMax, yMax) / 2.0;
long double rMin = 0.2 * min(xMax, yMax) / 2.0;
long double fanTurns = 2.0;
long double fanAngle = 0.0;
int printPercentageSteps = 20;
long double percentageStop = 100.0;

long double pi = atan(1) * 4;
long double F = 10.0;
long double fanWidth = 0.1;
bool printWork = false;
long double C_d = 2.0;
long double fanArea = dx * dy;
//todo: area should depende upon intersection
// of the element with the fan.
unsigned int stepsUntilPrint = 2000;
bool printPercentage = true;
bool onlyPrintFan = true;

