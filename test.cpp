#include <cmath>
#include <climits>
#include <Eigen/Dense>
#include "datagen.h"
#include "gtest/gtest.h"
#include "p35p_solver.h"
#include "p35p_solver_terminate.h"
#include "p35p_solver_initialize.h"
using namespace Eigen;

TEST(TestCaseName, TestName) {
	double eps = 1e-5;
	int it_num = 10; //iterations
	int succ_num = 0;
	double succ_threshold = 1e-5; //???
	double pass_threshold = 0.5;
	
	//allocate for generated data
	double* X = new double[12];
	double* x = new double[4];
	double* y = new double[4]; 
	double f_gen; Matrix3d R_gen; Vector3d C_gen;
	
	//allocate for estimated data
	double solution_num;
	double f_sol_data[10];
	int f_sol_size[2];
	double R_sol_data[90];
	int R_sol_size[3];
	double T_sol_data[30];
	int T_sol_size[2];
	p35p_solver_initialize();
	
	
	for(int curr_it = 0; curr_it < it_num; ++curr_it) {
		generateData(X, x, y, f_gen, R_gen, C_gen);

		p35p_solver(X, x, y, eps, &solution_num, f_sol_data,
				  f_sol_size, R_sol_data, R_sol_size, T_sol_data, T_sol_size);
		
		//allocate for comparison
		double min_diff = INT_MAX; //todo:ok??????????
		double diff_C = INT_MAX;
		double diff_R = INT_MAX;
		

		for (int i = 0; i < solution_num; ++i) { //todo: make some kind of round
			double f = f_sol_data[i];
			Vector3d T;
			Matrix3d R;
			for(int j = 0; j < 3; ++j) {
				T(j) = T_sol_data[3 * int(solution_num) + j];
			}
			int ind = 9*solution_num;
			for(int k = 0; k < 3; ++k) {
				for(int j = 0; j < 3; ++j) {
					R(j, k) = R_sol_data[ind]; //todo:check;
					ind++;
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
			} else if (diff == min_diff && diffR < diff_R) {
					diff_R = diffR;
					diff_C = diffC;
			}
		}
		if(min_diff < succ_threshold)
			succ_num++;
	}
	p35p_solver_terminate();
	ASSERT_GT(succ_num / it_num, pass_threshold);
	ASSERT_EQ(2, 3);
	delete[] X, x, y;
}