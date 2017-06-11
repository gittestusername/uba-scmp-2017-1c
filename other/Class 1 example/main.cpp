#include<iostream>
#include<vector>

using namespace std;

int main() {

    //Construct system:

    vector<double> curr;
    vector<double> last;


    int n = 400;
    int steps = 100;
    double r = 0.1;

    for (int i = 0; i < n; ++i) {
        if (i == n - 1) {
            last.push_back(0);
            curr.push_back(0);
        } else {
            last.push_back(i);
            curr.push_back(i);
        }

    }

    //Simulation:

    for (int i = 0; i < steps; ++i) {
        for (int i = 1; i < n - 1; ++i) {
            curr[i] = (1.0 - 2.0 * r) * last[i] + r * (last[i - 1] + last[i + 1]) / 2;
        }


        for (int i = 1; i < n - 1; ++i) {
            last[i] = curr[i];
        }


        //Print:
        cout << "[";
        for (int i = 0; i < n; ++i) {
            cout << curr[i] << ",";
        }
        cout << "],";
    }



}