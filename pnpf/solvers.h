//
// Created by elizaveta on 28.10.2019.
//

#ifndef PNP_TEST_SOLVERS_H
#define PNP_TEST_SOLVERS_H

#endif // PNP_TEST_SOLVERS_H

#include "pnpf.h"
#include <Eigen/Dense>
#include <climits>

using namespace Eigen;

template <class T> class Solver {
protected:
  T e;
  virtual void vSolve(const Matrix<T, 3, 4> &points_3d,
                      const Matrix<T, 2, 4> &points_2d, int *n, T *fs,
                      Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs, T diag) = 0;

public:
  Solver() { e = std::numeric_limits<T>::epsilon(); }
  void solve(const Matrix<T, 3, 4> &points_3d, const Matrix<T, 2, 4> &points_2d,
             int *n, T *fs, Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs,
             T diag = 1) {
    vSolve(points_3d, points_2d, n, fs, Rs, Cs, diag);
  };
};

template <class T> class MatlabSolver : public Solver<T> {
protected:
  T X[12];
  T x[4], y[4];
  int sol_num;
  int f_size[2], r_size[3], t_size[2];
  T f_data[10], r_data[90], t_data[30]; // todo: check sizes

  void dataEigenToMatlab(const Matrix<T, 3, 4> &points_3d,
                         const Matrix<T, 2, 4> &points_2d, T diag) {
    int ind = 0;
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
        X[ind] = points_3d(i, j);
        ind++;
      }
    }
    for (int i = 0; i < 4; ++i) {
      x[i] = points_2d(0, i) / diag;
      y[i] = points_2d(1, i) / diag;
    }
  }

  void dataMatlabToEigen(int *n, T *fs, Matrix<T, 3, 3> *Rs,
                         Matrix<T, 3, 1> *Cs, T diag) const {
    *n = sol_num;
    for (int i = 0; i < *n; ++i) {
      fs[i] = f_data[i] * diag;
      for (int j = 0; j < 3; ++j) {
        Cs[i](j) = t_data[3 * i + j];
      }
      int ind = 9 * i;
      for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
          Rs[i](j, k) = r_data[ind];
          ind++;
        }
      }
    }
  }
  void vSolve(const Matrix<T, 3, 4> &points_3d,
              const Matrix<T, 2, 4> &points_2d, int *n, T *fs,
              Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs, T diag) override {
    dataEigenToMatlab(points_3d, points_2d, diag);
    vMatlabSolve();
    dataMatlabToEigen(n, fs, Rs, Cs, diag);
  }

  virtual void vMatlabSolve() = 0;
};

template <class T> class P35PSolver : public MatlabSolver<T> {
protected:
  void vMatlabSolve(){};
};

template <> class P35PSolver<float> : public MatlabSolver<float> {
protected:
  void vMatlabSolve() override {
    p35pf_single(X, x, y, e, &sol_num, f_data, f_size, r_data, r_size, t_data,
                t_size);
  };
};

template <> class P35PSolver<double> : public MatlabSolver<double> {
protected:
  void vMatlabSolve() override {
    p35pf_double(X, x, y, e, &sol_num, f_data, f_size, r_data, r_size, t_data,
                t_size);
  };
};

template <class T> class P4PSolver : public MatlabSolver<T> {
protected:
  void vMatlabSolve(){};
};

template <> class P4PSolver<float> : public MatlabSolver<float> {
protected:
  void vMatlabSolve() override {
    p4pf_single(X, x, y, e, &sol_num, f_data, f_size, r_data, r_size, t_data,
                t_size);
  };
};

template <> class P4PSolver<double> : public MatlabSolver<double> {
protected:
  void vMatlabSolve() override {
    p4pf_double(X, x, y, e, &sol_num, f_data, f_size, r_data, r_size, t_data,
                t_size);
  };
};