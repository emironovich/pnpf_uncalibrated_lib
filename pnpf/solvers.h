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
                      Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs,
                      T diag) const = 0;

public:
  Solver() { e = std::numeric_limits<T>::epsilon(); }
  void solve(const Matrix<T, 3, 4> &points_3d, const Matrix<T, 2, 4> &points_2d,
             int *n, T *fs, Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs,
             T diag = 1) const {
    vSolve(points_3d, points_2d, n, fs, Rs, Cs, diag);
  };
};

template <class T> class P35PSolver : public Solver<T> {
protected:
  void vSolve(const Matrix<T, 3, 4> &points_3d,
              const Matrix<T, 2, 4> &points_2d, int *n, T *fs,
              Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs, T diag) const {};
};

template <> class P35PSolver<float> : public Solver<float> {
protected:
  void vSolve(const Matrix<float, 3, 4> &points_3d,
              const Matrix<float, 2, 4> &points_2d, int *n, float *fs,
              Matrix<float, 3, 3> *Rs, Matrix<float, 3, 1> *Cs,
              float diag) const override {

    float X[12];
    int ind = 0;
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
        X[ind] = points_3d(i, j);
        ind++;
      }
    }
    float x[4], y[4];
    for (int i = 0; i < 4; ++i) {
      x[i] = points_2d(0, i) / diag;
      y[i] = points_2d(1, i) / diag;
    }

    int f_size[2], r_size[3], t_size[2];
    float f_data[10], r_data[90], t_data[30]; // todo: check sizes

    p35p_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

    for (int i = 0; i < *n; ++i) {
      fs[i] = f_data[i] * diag;
      for (int j = 0; j < 3; ++j) {
        Cs[i](j) = t_data[3 * i + j];
      }
      ind = 9 * i;
      for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
          Rs[i](j, k) = r_data[ind];
          ind++;
        }
      }
    }
  };
};

template <> class P35PSolver<double> : public Solver<double> {
protected:
  void vSolve(const Matrix<double, 3, 4> &points_3d,
              const Matrix<double, 2, 4> &points_2d, int *n, double *fs,
              Matrix<double, 3, 3> *Rs, Matrix<double, 3, 1> *Cs,
              double diag) const override {

    double X[12];
    int ind = 0;
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
        X[ind] = points_3d(i, j);
        ind++;
      }
    }
    double x[4], y[4];
    for (int i = 0; i < 4; ++i) {
      x[i] = points_2d(0, i) / diag;
      y[i] = points_2d(1, i) / diag;
    }

    int f_size[2], r_size[3], t_size[2];
    double f_data[10], r_data[90], t_data[30]; // todo: check sizes

    p35p_double(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

    for (int i = 0; i < *n; ++i) {
      fs[i] = f_data[i] * diag;
      for (int j = 0; j < 3; ++j) {
        Cs[i](j) = t_data[3 * i + j];
      }
      ind = 9 * i;
      for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
          Rs[i](j, k) = r_data[ind];
          ind++;
        }
      }
    }
  };
};

template <class T> class P4PSolver : public Solver<T> {
protected:
  void vSolve(const Matrix<T, 3, 4> &points_3d,
              const Matrix<T, 2, 4> &points_2d, int *n, T *fs,
              Matrix<T, 3, 3> *Rs, Matrix<T, 3, 1> *Cs, T diag) const {};
};

template <> class P4PSolver<float> : public Solver<float> {
protected:
  void vSolve(const Matrix<float, 3, 4> &points_3d,
              const Matrix<float, 2, 4> &points_2d, int *n, float *fs,
              Matrix<float, 3, 3> *Rs, Matrix<float, 3, 1> *Cs,
              float diag) const override {
    float X[12];
    int ind = 0;
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
        X[ind] = points_3d(i, j);
        ind++;
      }
    }
    float x[4], y[4];
    for (int i = 0; i < 4; ++i) {
      x[i] = points_2d(0, i) / diag;
      y[i] = points_2d(1, i) / diag;
    }

    int f_size[2], r_size[3], t_size[2];
    float f_data[10], r_data[90], t_data[30]; // todo: check sizes

    p4pf_single(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

    for (int i = 0; i < *n; ++i) {
      fs[i] = f_data[i] * diag;
      for (int j = 0; j < 3; ++j) {
        Cs[i](j) = t_data[3 * i + j];
      }
      ind = 9 * i;
      for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
          Rs[i](j, k) = r_data[ind];
          ind++;
        }
      }
    }
  };
};

template <> class P4PSolver<double> : public Solver<double> {
protected:
  void vSolve(const Matrix<double, 3, 4> &points_3d,
              const Matrix<double, 2, 4> &points_2d, int *n, double *fs,
              Matrix<double, 3, 3> *Rs, Matrix<double, 3, 1> *Cs,
              double diag) const override {
    double X[12];
    int ind = 0;
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
        X[ind] = points_3d(i, j);
        ind++;
      }
    }
    double x[4], y[4];
    for (int i = 0; i < 4; ++i) {
      x[i] = points_2d(0, i) / diag;
      y[i] = points_2d(1, i) / diag;
    }

    int f_size[2], r_size[3], t_size[2];
    double f_data[10], r_data[90], t_data[30]; // todo: check sizes

    p4pf_double(X, x, y, e, n, f_data, f_size, r_data, r_size, t_data, t_size);

    for (int i = 0; i < *n; ++i) {
      fs[i] = f_data[i] * diag;
      for (int j = 0; j < 3; ++j) {
        Cs[i](j) = t_data[3 * i + j];
      }
      ind = 9 * i;
      for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
          Rs[i](j, k) = r_data[ind];
          ind++;
        }
      }
    }
  };
};