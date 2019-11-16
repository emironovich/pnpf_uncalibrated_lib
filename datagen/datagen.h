#pragma once
#include <Eigen/Dense>
#include <Eigen/LU>
#include <random>
#include <unsupported/Eigen/MatrixFunctions>

template <class Type>
Eigen::Matrix<Type, 3, 3> makeSkew(const Eigen::Matrix<Type, 3, 1> &a) {
  Eigen::Matrix<Type, 3, 3> S;
  S.setZero();
  S(0, 1) = -a(2);
  S(0, 2) = a(1);
  S(1, 2) = -a(0);

  S(1, 0) = a(2);
  S(2, 0) = -a(1);
  S(2, 1) = a(0);
  return S;
}

template <class Type>
void generateData(Eigen::Matrix<Type, 3, 4> &points_3d,
                  Eigen::Matrix<Type, 2, 4> &points_2d, Type &f,
                  Eigen::Matrix<Type, 3, 3> &R, Eigen::Matrix<Type, 3, 1> &C,
                  Type d = 0) {

  std::random_device dev;
  std::mt19937_64 generator(dev()); // Mersenne Twister generator
  std::uniform_real_distribution<Type> uniformDistribution(-1., 1.);
  auto uniform = [&]() { return uniformDistribution(generator); };

  // focal distance
  f = 200 + 1800 * (uniformDistribution(generator) + 1) / 2;

  // rotation
  Eigen::Matrix<Type, 3, 1> rVec =
      Eigen::Matrix<Type, 3, 1>::NullaryExpr(3, 1, uniform);
  rVec /= 2; //???
  Eigen::Matrix<Type, 3, 3> rVecSkew = makeSkew(rVec);
  R = rVecSkew.exp();

  // camera center
  C = Eigen::Matrix<Type, 3, 1>::NullaryExpr(3, 1, uniform);
  C /= 2; //???

  // calibration Eigen::Matrix
  Eigen::Matrix<Type, 3, 3> K;
  K.setZero();
  K(0, 0) = K(1, 1) = f;
  K(2, 2) = 1;

  // projection Eigen::Matrix
  Eigen::Matrix<Type, 3, 4> P;
  P.setIdentity();
  P.col(3) = -C;
  P = K * R * P;

  // points in space
  points_3d.setZero();
  points_3d.row(2).setConstant(6);
  Eigen::Matrix<Type, 3, 4> tmp =
      Eigen::Matrix<Type, 3, 4>::NullaryExpr(3, 4, uniform);
  points_3d += 2 * tmp;
  for (int i = 0; i < 4; ++i) {
    points_3d.col(i) = (R.transpose() * points_3d.col(i) + C).eval();
  }

  // image points
  Eigen::Matrix<Type, 4, 4> XMHom;
  XMHom << points_3d, 1, 1, 1, 1;
  Eigen::Matrix<Type, 3, 4> pHom = P * XMHom;

  for (int i = 0; i < 4; ++i) {
    Eigen::Matrix<Type, 2, 1> dist;
    if (d > 0) {
      std::normal_distribution<Type> normalDistribution(0, d);
      auto normal = [&]() { return normalDistribution(generator); };
      dist = Eigen::Matrix<Type, 2, 1>::NullaryExpr(2, 1, normal);
    } else {
      dist.setZero();
    }
    points_2d.col(i) = pHom.col(i).hnormalized() + dist;
  }
}