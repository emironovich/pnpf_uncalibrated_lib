#pragma once
#include <Eigen/Dense>
using namespace Eigen;

Matrix3d makeSkew(Vector3d& a);
void generateData(double** X, double* x, double* y, double& f, Matrix3d& R, Vector3d& C, double d = 0);
