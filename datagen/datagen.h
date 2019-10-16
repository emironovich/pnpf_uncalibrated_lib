#pragma once
#include <Eigen/Dense>
#include <Eigen/LU>
#include <random>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;

template <class Type> Matrix<Type, 3, 3> makeSkew(const Matrix<Type, 3, 1> &a) {
  Matrix<Type, 3, 3> S;
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
void generateData(Type *X, Type *x, Type *y, Type &f, Matrix<Type, 3, 3> &R,
                  Matrix<Type, 3, 1> &C, Type d = 0) {
  // focal distance
  f = 200 + 1800 * (Type)(rand()) / RAND_MAX; // todo: change function

  // rotation
  Matrix<Type, 3, 1> rVec;
  rVec.setRandom();
  rVec /= 2; //???
  Matrix<Type, 3, 3> rVecSkew = makeSkew(rVec);
  R = rVecSkew.exp();

  // camera center
  C.setRandom();
  C /= 2; //???

  // calibration matrix
  Matrix<Type, 3, 3> K;
  K.setZero();
  K(0, 0) = K(1, 1) = f;
  K(2, 2) = 1;

  // projection matrix
  Matrix<Type, 3, 4> P;
  P.setIdentity();
  P.col(3) = -C;
  P = K * R * P;

  // points in space
  Matrix<Type, 3, 4> XM;
  XM.setZero();
  XM.row(2).setConstant(6);
  Matrix<Type, 3, 4> tmp;
  tmp.setRandom();
  XM += 2 * tmp;
  for (int i = 0; i < 4; ++i) {
    XM.col(i) = (R.transpose() * XM.col(i) + C).eval();
  }
  int ind = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 3; ++i) { // assuming matlab made column-major array
      X[ind] = XM(i, j);
      ind++;
    }
  }

  /*for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 4; ++j) {
          X[i][j] = XM(i, j);
      }
  }*/

  // image points
  Matrix<Type, 4, 4> XMHom;
  XMHom << XM, 1, 1, 1, 1;
  Matrix<Type, 3, 4> pHom = P * XMHom;

  for (int i = 0; i < 4; ++i) {
    Type xDist, yDist;
    if (d > 0) {
      std::default_random_engine generator;
      std::normal_distribution<Type> distribution(0, d);
      xDist = distribution(generator);
      yDist = distribution(generator);
    } else {
      xDist = yDist = 0;
    }
    x[i] = pHom(0, i) / pHom(2, i) + xDist;
    y[i] = pHom(1, i) / pHom(2, i) + yDist;
  }
}