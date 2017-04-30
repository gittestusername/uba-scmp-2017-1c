#ifndef MAT_H
#define MAT_H
#include <iostream>
#include <vector>
#include <math.h>
#include <sstream>
#include <cmath>   

long double pi = acos(-1.0L);

template <class T> class mat {
  private:
    std::vector<T> data_;
    int cols_;
    int rows_;

  public:
    mat(int rows, int cols);
    mat(int rows, int cols, T value);
    ~mat();
    T at(int rows, int cols);
    T at(int rows);
    void set(int rows, int cols, T value);
    void set(int rows, T value);

    void print();
    mat<T> prod(mat m);
    int cols();
    int rows();
    void addAt(int rows, int cols, T value);
    void id();
    mat<T> operator+(mat<T>&);
    void operator+=(mat<T>&);
    mat<T> operator-(mat<T>&);
    void operator-=(mat<T>&);
    mat<T> operator*(mat<T>&);
    mat<T> jacobi(mat<T> b);
    mat<T> gaussElimination(mat<T> b);
    mat<T> copy();



};

template<class T> inline std::ostream& operator<<(std::ostream& os, mat<T>& m) {
    for (int i = 0; i < m.rows(); ++i) {
        std::cout << std::endl;
        for (int j = 0; j < m.cols(); ++j) {
            std::stringstream ss;
            ss << m.at(i, j);
            std::string str = ss.str();
            int len = str.length();
            int maxSpaces = 7; //TODO: Use max of matrix instead of a fixed value.
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

template <class T> void mat<T>::print() {
    for (int i = 0; i < rows_; ++i) {
        std::cout << std::endl;
        for (int j = 0; j < cols_; ++j) {
            std::stringstream ss;
            ss << data_[i * cols_ + j];
            std::string str = ss.str();
            int len = str.length();
            int maxSpaces = 7; //TODO: Use max of matrix instead of a fixed value.
            int spaces = maxSpaces - len;
            while (spaces > 0) {
                std::cout << " ";
                spaces--;
            }
            std::cout << " " << data_[i * cols_ + j];
        }
    }
    std::cout << std::endl;
}


template <class T> T mat<T>::at(int rows, int cols) {
    return data_[rows * cols_ + cols];
}


template <class T> T mat<T>::at(int rows) {
    return data_[rows];
}

template <class T> void mat<T>::set(int rows, int cols, T value) {
    data_[rows * cols_ + cols] = value;
}


template <class T> void mat<T>::set(int rows, T value) {
    data_[rows] = value;
}


template<typename T>
mat<T> mat<T>::operator*( mat<T>& m) {
    int mr = m.rows();
    int mc = m.cols();
    if (cols_ != mr) {
        //ERROR. //TODO: Implement.
    }
    mat<T> res(rows_, mc, 0);

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < mc; ++j) {
            for (int k = 0; k < rows_; ++k) {
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



template <class T> mat<T> mat<T>::jacobi(mat<T> b) {
//TODO: Custom seed, threshold and niters version overloading.
//TODO: Catch size problems function.
    mat<T> x(rows(), 1, 0);
    int nIters = 30;
    for (int itr = 0; itr < nIters; ++itr)
    {
        for (int i = 0; i < cols(); ++i)
        {
            T sum;
            sum = 0;
            for (int j = 0; j < cols(); ++j)
            {
                if(i == j) continue;
                sum += at(i,j)*x.at(j);
            }
            x.set(i, (b.at(i)-sum)/at(i,i));
        }
    }
    return x;
}


template <class T> mat<T> mat<T>::copy() {
    mat<T> res(rows(), cols(), 0);
    for (int i = 0; i < rows(); ++i)
    {
        for (int j = 0; j < cols(); ++j)
        {
            res.set(i,j,at(i,j));
        }
    }
    return res;
}





template <class T> mat<T> mat<T>::gaussElimination(mat<T> b) {
//TODO: Custom seed, threshold and niters version overloading.
//TODO: Catch size problems function.
//Warning: SLOW.
    mat<T> X(rows(), 1, 0);
    mat<T> augA(rows(),cols()+1, 0);
    augA += *this;
    for (int i = 0; i < augA.rows(); ++i)
    {
        augA.set(i, augA.cols()-1, b.at(i));
    }


    for (int i=0; i< cols(); i++) {
        // Search for maximum in this column
        T maxElem = std::abs(augA.at(i,i));
        int maxRow = i;
        for (int k = i+1; k < cols(); k++) {
            if (std::abs(augA.at(k,i) > maxElem)) {
                maxElem = std::abs(augA.at(k,i));
                maxRow = k;
            }
        }
        if(maxElem == 0){
            std::cout << "ERROR: Matrix is singular" << std::endl;
            mat<T> er(rows(), 0);
            return er;
            //TODO: Using the zero vector as a
            // result and printing ERROR, tries
            // to make the mistake visible,
            // but we should, instead, rise
            // an exeption.
        }

 

        // Swap maximum row with current row (column by column)
        for (int k=i; k<rows()+1;k++) {
            T tmp = augA.at(maxRow,k);
            augA.set(maxRow, k, augA.at(i,k));
            augA.set(i, k, tmp);
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<rows(); k++) {
            T c = -augA.at(k, i)/augA.at(i,i);
            for(int j=i; j<rows()+1; j++) {
                if (i==j) {
                    augA.set(k,j,0);
                } else {
                    augA.addAt(k, j, c*augA.at(i,j));
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i=cols()-1; i >= 0; i--) {
        X.set(i, augA.at(i,cols())/augA.at(i,i));
        for (int k = i-1; k >= 0; k--) {
            augA.set(k,cols(), augA.at(k,cols()) - augA.at(k,i) * X.at(i));
        }
    }
    return X;
}








#endif


