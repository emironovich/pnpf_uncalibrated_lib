/*
#ifndef PNP_TEST_TEST_H
#define PNP_TEST_TEST_H

#endif //PNP_TEST_TEST*/

#include <cmath>
#include <climits>
#include <Eigen/Dense>
#include "datagen.h"
#include "p35p_solver.h"
#include "p4p_solver.h"

using namespace Eigen;

struct TestResult{
    int existSolutions = 0;
    int belowThreshold = 0;
};

template <class Type>
TestResult runFunction(const char * funcName, Type eps, int it_num = 1) {
    //Type eps = 1e-5;
    //int it_num = 10000; //iterations
    int succ_num = 0;
    Type succ_threshold = 1.5e-9; //todo: decide if want to be an argument
    Type pass_threshold = 0.9;
    int zero_solutions_num = 0;

    //allocate for generated data
    Type X[12];
    Type x[4];
    Type y[4];
    Type f_gen; Matrix<Type, 3, 3> R_gen; Matrix<Type, 3, 1> C_gen;

    //allocate for estimated data
    Type solution_num = 0;
    Type f_sol_data[10];
    int f_sol_size[2];
    Type R_sol_data[90];
    int R_sol_size[3];
    Type T_sol_data[30];
    int T_sol_size[2];
    //p4p_solver_initialize();


    for(int curr_it = 0; curr_it < it_num; ++curr_it) {
        generateData(X, x, y, f_gen, R_gen, C_gen);

        if(strcmp(funcName, "p4p"))
            p4p_solver(X, x, y, eps, &solution_num, f_sol_data,
                                      f_sol_size, R_sol_data, R_sol_size, T_sol_data, T_sol_size);
        else if(strcmp(funcName, "p3.5p"))
            p35p_solver(X, x, y, eps, &solution_num, f_sol_data,
                                          f_sol_size, R_sol_data, R_sol_size, T_sol_data, T_sol_size);
        else {
            solution_num = 0; //todo: maybe add a warning
        }

        //allocate for comparison
        Type min_diff = INT_MAX; //todo:ok??????????
        Type diff_C = INT_MAX;
        Type diff_R = INT_MAX;


        for (int i = 0; i < solution_num; ++i) { //todo: make some kind of round
            Type f = f_sol_data[i];
            Matrix<Type, 3, 1>  T;
            Matrix<Type, 3, 3> R;
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

            Matrix<Type, 3, 3> K;
            K.setZero();
            K(0, 0) = K(1, 1) = f;
            K(2, 2) = 1;

            Matrix<Type, 3, 1>  C = (-R).colPivHouseholderQr().solve(T); //todo: fix

            Type diff = abs(f - f_gen);

            Type diffR = (R-R_gen).norm()/3;
            Type diffC = (C - C_gen).norm()/C_gen.norm();
            if (diff < min_diff) {
                min_diff = diff;
                diff_R = diffR;
                diff_C = diffC;
            } else if (diff == min_diff && diffR < diff_R) {
                diff_R = diffR;
                diff_C = diffC;
            }
        }
        if(min_diff / f_gen < succ_threshold && solution_num != 0)
            succ_num++;
        if(solution_num == 0)
            zero_solutions_num++;
    }

    if(strcmp(funcName, "p4p"))
        p4p_solver_terminate();
    else if(strcmp(funcName, "p3.5p"))
        p35p_solver_terminate();

    return {(it_num - zero_solutions_num) , succ_num};
    //ASSERT_GT(Type(succ_num) / (it_num - zero_solutions_num) , pass_threshold);
}