#include "pch.h"
#include <cmath>
#include <Eigen/Dense>
#include "datagen.h"
using namespace Eigen;

TEST(TestCaseName, TestName) {
	double eps = 1e-5;
	
	double** X = new double*[3];
	for (int i = 0; i < 3; ++i) {
		X[i] = new double[4];
	}
	double* x = new double[4];
	double* y = new double[4]; 
	double f_gen; Matrix3d R_gen; Vector3d C_gen;


	generateData(X, x, y, f_gen, R_gen, C_gen);
	
	//apply function
	//extern void p35p_solver(const float X[12], const float x[4], const float y[4],
  //float e, float *solution_num, float f_sol_data[], int f_sol_size[2], float
 // R_sol_data[], int R_sol_size[3], float T_sol_data[], int T_sol_size[2]);
 
	int solution_num;
	
	for (int i = 0; i < solution_num; ++i) {
		double f = f_sol[i];
		Vector3d T;
		Matrix3d R;
		for(int j = 0; j < 3; ++j) {
			T(j) = T_sol[j][i];
			for(int k = 0; k < 3; ++k){
				R(j, k) = R_sol(j, k, i); //check dimentions
			}
		}
		Matrix3d K;
		K.setZero();
		K(0, 0) = K(1, 1) = f;
		K(2, 2) = 1;
		
		Vector3d C = (-K*R).colPivHouseholderQr().solve(T);
		
		double diff = abs(f - f_gen);

        double diffR = (R-R_gen).norm()/3;
        double diffC = (C - C_gen).norm()/C_gen.norm();
		if (diff < min_diff) {
                min_diff = diff;    
                diff_R = diffR;
                diff_C = diffC;
		} elseif (diff == min_diff && diffR < diff_R) {
                diff_R = diffR;
                diff_C = diffC;
        }	
	}
	//assert?
	
	for (int i = 0; i < 3; ++i) {
		delete[] X[i];
	}
	delete[] X, x, y;
}