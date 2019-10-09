#pragma once
#include <Eigen/Dense>
using namespace Eigen;

Matrix3d makeSkew(Vector3d& a);
Matrix3f makeSkew(Vector3f& a);
void generateData(double* X, double* x, double* y, double& f, Matrix3d& R, Vector3d& C, double d = 0);
void generateData(float* X, float* x, float* y, float& f, Matrix3f& R, Vector3f& C, float d = 0);