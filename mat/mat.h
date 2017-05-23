#ifndef MAT_H
#define MAT_H
#include <iostream>
#include <vector>
#include <math.h>
#include <sstream>
#include <cmath>

using std::cout;
using std::endl;
long double pi = acos(-1.0L);

template <class T> class mat {
  private:
    std::vector<T> data_;
    int cols_;
    int rows_;
    void addRowTimes(int a, int b, T k);
    void divRowTimes(int row, T k);

    
  public:
    mat(int rows, int cols);
    mat(int rows, int cols, T value);
    ~mat();
    T at(int rows, int cols);
    T at(int rows);
    void set(int rows, int cols, T value);
    void set(int rows, T value);

    int cols();
    int rows();
    void addAt(int rows, int cols, T value);
    void addAt(int rows, T value);
    void setRow(int row, T value);
    void setCol(int col, T value);

    void id();
    mat<T> operator+(mat<T>&);
    void operator+=(mat<T>&);
    mat<T> operator-(mat<T>&);
    void operator-=(mat<T>&);
    mat<T> operator*(mat<T>&);

    mat<T> jacobi(mat<T> b);
    mat<T> sparseProd(mat<T> m);
    mat<T> gaussElimination(mat<T> b);
    mat<T> copy();
    void fillScheme(mat<T> scheme);
    mat<T> inv();

//Para el metodo fuertemente implicito, el de dos pasos, se resuelve asi: uj+1 = jacobi(A,(C*uj))

};

template<class T> inline std::ostream& operator<<(std::ostream& os, mat<T>& m) {
    for (int i = 0; i < m.rows(); ++i) {
        std::cout << std::endl;
        for (int j = 0; j < m.cols(); ++j) {
            std::stringstream ss;
            ss << m.at(i, j);
            std::string str = ss.str();
            int len = str.length();
            int maxSpaces = 10; //TODO: Use max of matrix instead of a fixed value.
            int spaces = maxSpaces - len;
            while (spaces > 0) {
                os << " ";
                spaces--;
            }
            os << " " << m.at(i, j);
        }
    }
    os << std::endl;
}

template<typename T>
mat<T> mat<T>::operator+( mat<T>& m) {
    mat<T> res(rows_, cols_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            res.set(i, j, at(i, j) + m.at(i, j));
        }
    }
    return res;
}

template<typename T>
void mat<T>::operator+=( mat<T>& m) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            set(i, j, at(i, j) + m.at(i, j));
        }
    }
}

template<typename T>
mat<T> mat<T>::operator-( mat<T>& m) {
    mat<T> res(rows_, cols_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            res.set(i, j, at(i, j) - m.at(i, j));
        }
    }
    return res;
}

template<typename T>
void mat<T>::operator-=( mat<T>& m) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            set(i, j, at(i, j) - m.at(i, j));
        }
    }
}

template <class T> mat<T>::mat(int rows, int cols) {
    std::vector<T> tmp(cols * rows);
    data_ = tmp; //TODO: Does this take O(nm) time?
    rows_ = rows;
    cols_ = cols;
}

template <class T> mat<T>::mat(int rows, int cols, T value) {
    std::vector<T> tmp(cols * rows);
    data_ = tmp;
    rows_ = rows;
    cols_ = cols;
    for (int i = 0; i < rows_ * cols_; ++i) {
        data_[i] = value;
    }
}

template <class T> mat<T>::~mat() {
//delete(data_); //Throws double free, but why?
}

template <class T> void mat<T>::id() {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            if (i == j) {
                set(i, j, 1);
            } else {
                set(i, j, 0);
            }
        }

    }
}




template <class T> void mat<T>::setCol(int col, T value) {
    for (int i = 0; i < rows(); ++i) {
        set(i, col, value);
    }
}




template <class T> void mat<T>::setRow(int row, T value) {
    for (int i = 0; i < cols(); ++i) {
        set(row, i, value);
    }
}


template <class T> T mat<T>::at(int row, int col) {
    return data_[row * cols_ + col];
}


template <class T> T mat<T>::at(int row) {
    if (cols() > 1 && this->rows() > 1)
        std::cout << " ERROR: Unidimentional indexing in multidimentional matrix" << std::endl;
    if (cols() <= row && this->rows() <= row){
        std::cout << " ERROR: Index out of Matrix" << std::endl;
        return 0;
    }
    return data_[row]; //Because its unidimentional, it doesnt matter if its vertical or horizontal.
}

template <class T> void mat<T>::set(int row, int col, T value) {
    if(row >= rows() || col >= cols())
        std::cout << " ERROR: Index out of Matrix" << std::endl;
    data_[row * cols() + col] = value;
}


template <class T> void mat<T>::set(int row, T value) {
    if (cols() > 1 && this->rows() > 1)
        std::cout << " ERROR: Unidimentional indexing in multidimentional matrix" << std::endl;
    data_[row] = value; //Because its unidimentional, it doesnt matter if its vertical or horizontal.
}


template<typename T>
mat<T> mat<T>::operator*( mat<T>& m) {
    int mr = m.rows();
    int mc = m.cols();
    if (cols() != mr) {
        //ERROR. //TODO: Implement.
    }
    mat<T> res(rows(), mc, 0);

    for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < mc; ++j) {
            for (int k = 0; k < cols(); ++k) {
                res.addAt(i, j, at(i, k)*m.at(k, j));
            }
        }
    }

    return res;
}

template <class T> int mat<T>::cols() {
    return cols_;
}

template <class T> int mat<T>::rows() {
    return rows_;
}

template <class T> void mat<T>::addAt(int rows, int cols, T value) {
    data_[rows * cols_ + cols] += value;
}

template <class T> void mat<T>::addAt(int rows, T value) {
    data_[rows] += value;
}

template <class T> mat<T> mat<T>::jacobi(mat<T> b) {
//TODO: Custom seed, threshold and n-iters version overloading.
//TODO: Catch size problems function.
    mat<T> x(rows(), 1, 0);
    int nIters = 25;
    for (int itr = 0; itr < nIters; ++itr) {
        for (int i = 0; i < cols(); ++i) {
            T sum;
            sum = 0;
            for (int j = 0; j < cols(); ++j) {
                if (i == j) continue;
                sum += at(i, j) * x.at(j);
            }
            x.set(i, (b.at(i) - sum) / at(i, i));
        }
    }
    return x;
}


template <class T> mat<T> mat<T>::copy() {
    mat<T> res(rows(), cols(), 0);
    for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
            res.set(i, j, at(i, j));
        }
    }
    return res;
}





template <class T> mat<T> mat<T>::gaussElimination(mat<T> b) {
//TODO: Custom seed, threshold and niters version overloading.
//TODO: Catch size problems function.
//Warning: SLOW.
    mat<T> X(rows(), 1, 0);
    mat<T> augA(rows(), cols() + 1, 0);
    augA += *this;
    for (int i = 0; i < augA.rows(); ++i) {
        augA.set(i, augA.cols() - 1, b.at(i));
    }


    for (int i = 0; i < cols(); i++) {
        // Search for maximum in this column
        T maxElem = std::abs(augA.at(i, i));
        int maxRow = i;
        for (int k = i + 1; k < cols(); k++) {
            if (std::abs(augA.at(k, i) > maxElem)) {
                maxElem = std::abs(augA.at(k, i));
                maxRow = k;
            }
        }
        if (maxElem == 0) {
            std::cout << "ERROR: Matrix is singular" << std::endl;
            mat<T> er(rows(), 1);
            return er;
            //TODO: Using the zero vector as a
            // result and printing ERROR, tries
            // to make the mistake visible,
            // but we should, instead, rise
            // an exeption.
        }



        // Swap maximum row with current row (column by column)
        for (int k = i; k < rows() + 1; k++) {
            T tmp = augA.at(maxRow, k);
            augA.set(maxRow, k, augA.at(i, k));
            augA.set(i, k, tmp);
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < rows(); k++) {
            T c = -augA.at(k, i) / augA.at(i, i);
            for (int j = i; j < rows() + 1; j++) {
                if (i == j) {
                    augA.set(k, j, 0);
                } else {
                    augA.addAt(k, j, c * augA.at(i, j));
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i = cols() - 1; i >= 0; i--) {
        X.set(i, augA.at(i, cols()) / augA.at(i, i));
        for (int k = i - 1; k >= 0; k--) {
            augA.set(k, cols(), augA.at(k, cols()) - augA.at(k, i) * X.at(i));
        }
    }
    return X;
}


template <class T> void mat<T>::fillScheme(mat<T> scheme) {
    if (  (scheme.rows() > 1) || (scheme.cols() % 2 != 1) ) {
        std::cout << "ERROR: Scheme matrix must be a one dimensional matrix of odd length" << std::endl;
        return;
    }

    //First we set everything to zero.
    for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
            set(i, j, 0);
        }
    }

    //Then we must fill the discretization into the matrix.
    int lenTails = (scheme.cols() - 1) / 2; //'(123)' (4) (567) = '(tail)' (mid) (tail)

    for (int i = 0; i < rows() ; ++i) {
        for (int j = 0; j < scheme.cols(); ++j) {
            int c = i - lenTails + j;
            if (c >= 0 && c < cols()) {
                set(i, c, scheme.at(j));
            }
        }
    }

    //Note that the boundary conditions are not trated here.
    return;
}

template<typename T>
mat<T> mat<T>::sparseProd(mat<T> m) {
//Halves the time needed for calculation if less than 5%
// of the elements are not zero.
    if (cols_ != m.rows()) {
        //ERROR. //TODO: Implement.
    }
    mat<T> res(rows(), m.cols(), 0);

    std::vector< std::vector< T > > rws;
    for (int i = 0; i < rows(); i++) {
        std::vector< T > r;
        for (int j = 0; j < cols(); j++) {
            if (at(i, j) != 0)
                r.push_back(j);
        }
        rws.push_back(r);
    }

    for (int mc = 0; mc < m.cols(); ++mc) {
        for (int i = 0; i < rows(); ++i) {
            for (int jElem = 0; jElem < rws[i].size(); ++jElem) {
                res.addAt(i, mc, at(i, rws[i][jElem])*m.at(rws[i][jElem], mc));
            }
        }
    }

    return res;
}

#endif




template<typename T>
mat<T> mat<T>::inv() {
    //TODO: Too slow. Better method? Inplace maybe?

    mat<T> tmp = this->copy();

    mat<T> inv(rows(), cols());
    inv.id();


    if(cols() != rows()){
        std::cout << "Error: Trying to invert a non square matrix" << std::endl;
        //TODO: Better error system.
        return inv;
    }


    int n = cols();



    for (int d = 0; d < n; ++d)
    {
        for (int r = 0; r < n; ++r)
        {
            if(r == d) continue;
            inv.addRowTimes(r, d, -tmp.at(r,d)/tmp.at(d,d));
            tmp.addRowTimes(r, d, -tmp.at(r,d)/tmp.at(d,d));

        }
    }

        for (int d = 0; d < n; ++d)
    {
        T div = tmp.at(d,d);
        tmp.divRowTimes(d, div);
        inv.divRowTimes(d, div);
    }

    return inv;
}


template<typename T>
void mat<T>::addRowTimes(int a, int b, T k){
    for (int idx = 0; idx < cols(); ++idx)
    {
        set(a ,idx, at(a,idx) + at(b,idx)*k);
    }
}



template<typename T>
void mat<T>::divRowTimes(int row, T k){
    for (int idx = 0; idx < cols(); ++idx)
    {
        set(row, idx, at(row,idx)/k);
    }
}