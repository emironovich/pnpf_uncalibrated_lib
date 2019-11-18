//
// Created by elizaveta on 28.10.2019.
//

#ifndef PNPF_SOLVERS_HPP
#define PNPF_SOLVERS_HPP

#include <Eigen/Dense>
#include <climits>

#include "pnpf/pnpf.h"

namespace pnpf {

template <typename T> struct SolverTraits {};

template <class Solver> class MatlabSolver {
public:
  using Scalar = typename SolverTraits<Solver>::Scalar;
  static constexpr int MaxSolutions = SolverTraits<Solver>::MaxSolutions;
  static constexpr Scalar epsilon = std::numeric_limits<Scalar>::epsilon();
  using Matrix34 = Eigen::Matrix<Scalar, 3, 4>;
  using Matrix24 = Eigen::Matrix<Scalar, 2, 4>;
  using Matrix33 = Eigen::Matrix<Scalar, 3, 3>;
  using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

  void dataEigenToMatlab(const Matrix34 &points_3d, const Matrix24 &points_2d,
                         const Scalar diag) {
    X = points_3d.data();
    Eigen::Map<Matrix24> xy_map(xy);
    xy_map = points_2d.array() * Scalar(1. / diag);
  }

  void dataMatlabToEigen(int *n, Scalar *fs, Matrix33 *Rs, Vector3 *Cs,
                         Scalar diag) {
    *n = sol_num;
    for (int i = 0; i < *n; ++i) {
      fs[i] = f_data[i] * diag;
      Cs[i] = Eigen::Map<Vector3>(t_data + 3 * i);
      Rs[i] = Eigen::Map<Matrix33>(r_data + 9 * i);
    }
  }
  void solve(const Matrix34 &points_3d, const Matrix24 &points_2d, int *n,
             Scalar *fs, Matrix33 *Rs, Vector3 *Cs, Scalar diag = 1) {
    dataEigenToMatlab(points_3d, points_2d, diag);
    static_cast<Solver *>(this)->solveImpl();
    dataMatlabToEigen(n, fs, Rs, Cs, diag);
  }

protected:
  const Scalar *X;
  Scalar xy[8];
  int sol_num;
  int f_size[2], r_size[3], t_size[2];
  Scalar f_data[MaxSolutions], r_data[MaxSolutions * 9],
      t_data[MaxSolutions * 3];
};

template <class Scalar> class P35PSolver {};
template <typename T> struct SolverTraits<P35PSolver<T>> {
  using Scalar = T;
  static constexpr int MaxSolutions = 10;
};

template <> class P35PSolver<float> : public MatlabSolver<P35PSolver<float>> {
public:
  void solveImpl() {
    p35pf_single(X, xy, epsilon, &sol_num, f_data, f_size, r_data, r_size,
                 t_data, t_size);
  };
};

template <> class P35PSolver<double> : public MatlabSolver<P35PSolver<double>> {
public:
  void solveImpl() {
    p35pf_double(X, xy, epsilon, &sol_num, f_data, f_size, r_data, r_size,
                 t_data, t_size);
  };
};

template <class Scalar> class P4PSolver {};
template <typename T> struct SolverTraits<P4PSolver<T>> {
  using Scalar = T;
  static constexpr int MaxSolutions = 10;
};

template <> class P4PSolver<float> : public MatlabSolver<P4PSolver<float>> {
public:
  void solveImpl() {
    p4pf_single(X, xy, epsilon, &sol_num, f_data, f_size, r_data, r_size,
                t_data, t_size);
  };
};

template <> class P4PSolver<double> : public MatlabSolver<P4PSolver<double>> {
public:
  void solveImpl() {
    p4pf_double(X, xy, epsilon, &sol_num, f_data, f_size, r_data, r_size,
                t_data, t_size);
  };
};

} // namespace pnpf

#endif // PNPF_H
