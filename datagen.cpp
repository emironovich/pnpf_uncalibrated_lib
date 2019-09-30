#include <Eigen/Dense>
#include <Eigen/LU> 
#include <unsupported/Eigen/MatrixFunctions>
#include <random>

using namespace Eigen;

Matrix3d makeSkew(Vector3d& a) {
	Matrix3d S;
	S.setZero();
	S(0, 1) = -a(2);
	S(0, 2) = a(1);
	S(1, 2) = -a(0);

	S(1, 0) = a(2);
	S(2, 0) = -a(1);
	S(2, 1) = a(0);
	return S;
}

void generateData(double** X, double* x, double* y, double& f, Matrix3d& R, Vector3d& C, double d = 0) {
	//focal distance
	f = 200 + 1800 * (double)(rand()) / RAND_MAX;

	//rotation
	Vector3d rVec;
	rVec.setRandom();
	rVec /= 2; //???
	Matrix3d rVecSkew = makeSkew(rVec);
	R = rVecSkew.exp();

	//camera center
	C.setRandom();
	C /= 2; //???

	//calibration matrix
	Matrix3d K;
	K.setZero();
	K(0, 0) = K(1, 1) = f;
	K(2, 2) = 1;

	//projection matrix
	Matrix<double, 3, 4> P;
	P.setIdentity();
	P.col(3) = -C;
	P = K * R * P;

	//points in space
	//todo: add eval??
	Matrix<double, 3, 4> XM;
	XM.setZero();
	XM.row(2).setConstant(6);
	XM += 2 * (MatrixXd::Random(3, 4));
	for (int i = 0; i < 4; ++i) {
		XM.col(i) = (R.transpose() * XM.col(i) + C).eval();
	}
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 4; ++j) {
			X[i][j] = XM(i, j);
		}
	}

	//image points
	Matrix4d XMHom;
	XMHom << XM,
		MatrixXd::Constant(1, 4, 1);
	Matrix<double, 3, 4> pHom = P * XMHom;

	for (int i = 0; i < 4; ++i) {
		double xDist, yDist;
		if (d > 0) {
			std::default_random_engine generator;
			std::normal_distribution<double> distribution(0, d);
			xDist = distribution(generator);
			yDist = distribution(generator);
		}
		else {
			xDist = yDist = 0;
		}
		x[i] = pHom(0, i) / pHom(2, i) + xDist;
		y[i] = pHom(0, i) / pHom(2, i) + yDist;
	}

}

