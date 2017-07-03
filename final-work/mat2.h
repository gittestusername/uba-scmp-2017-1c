#include <iostream>
using namespace std;

class mat2 {
  public:
    long double * data;
    int iMax;
    int jMax;
    mat2() {
        iMax = 0;
        jMax = 0;
        data = NULL;
    }

    mat2(int iMax, int jMax) {
        this->iMax = iMax;
        this->jMax = jMax;
        inic(0);
    }


    mat2(int iMax, int jMax, long double val) {
        this->iMax = iMax;
        this->jMax = jMax;
        inic(val);
    }

    void inic(long double val) {
        int n = iMax * jMax;
        data = new long double[n];

        for (int k = 0; k < n; k++)
            data[k] = val;
    }

    void set(int i, int j, long double val) {
        data[i * jMax + j] = val;
    }

    void setRow(int i, mat2 &B) {
        for (int j = 0; j < cols(); ++j) {
            set(i, j, B.at(0, j));
        }
    }

    void setAll(mat2 &B) {
        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                set(i, j, B.at(i, j));
            }
        }
    }


    void add(int i, int j, long double val) {
        data[i * jMax + j] += val;
    }


    long double at(int i, int j) {
        return data[i * jMax + j];
    }


    long double cols() {
        return jMax;
    }


    long double rows() {
        return iMax;
    }


    void print() {
        for (int i = 0; i < iMax; i++) {
            for (int j = 0; j < jMax; j++) {
                cout << data[i * jMax + j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }





    ~mat2() {
        if (data != NULL)
            delete [] data;
    }
};