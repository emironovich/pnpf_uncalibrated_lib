#include <cmath>
#include <climits>
#include <Eigen/Dense>
#include "datagen.h"
#include "gtest/gtest.h"
#include "p35p_solver_double.h"
#include "p35p_solver_single.h"
#include "solve_P4Pf_double.h"
#include "solve_P4Pf_mixed.h"
//#include "p35p_solver_terminate.h"
//#include "p35p_solver_initialize.h"
using namespace Eigen;

TEST(AccuracyCheck, P35PDouble) {
    double eps = 1e-5;
    int it_num = 10000; //iterations
    int succ_num = 0;
    double succ_threshold = 1.7e-11;
    double pass_threshold = 0.9;
    int zero_solutions_num = 0;

    //allocate for generated data
    double X[12];
    double x[4];
    double y[4];
    double f_gen; Matrix3d R_gen; Vector3d C_gen;

    //allocate for estimated data
    double solution_num;
    double f_sol_data[10];
    int f_sol_size[2];
    double R_sol_data[90];
    int R_sol_size[3];
    double T_sol_data[30];
    int T_sol_size[2];
    //p35p_solver_initialize();


    for(int curr_it = 0; curr_it < it_num; ++curr_it) {
        generateData(X, x, y, f_gen, R_gen, C_gen);

        p35p_double::p35p_solver_double(X, x, y, eps, &solution_num, f_sol_data,
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
        if(min_diff / f_gen < succ_threshold && solution_num != 0)
            succ_num++;
        if(solution_num == 0)
            zero_solutions_num++;
    }
    p35p_double::p35p_solver_double_terminate();
    ASSERT_GT(double(succ_num) / (it_num - zero_solutions_num) , pass_threshold);
    //todo: add check for T and C
}


TEST(AccuracyCheck, P35PSingle) {
    float eps = 1e-5;
    int it_num = 10000; //iterations
    int succ_num = 0;
    float succ_threshold = 0.006;
    float pass_threshold = 0.9;
    int zero_solutions_num = 0;

    //allocate for generated data
    float X[12];
    float x[4];
    float y[4];
    float f_gen; Matrix3f R_gen; Vector3f C_gen;

    //allocate for estimated data
    float solution_num;
    float f_sol_data[10];
    int f_sol_size[2];
    float R_sol_data[90];
    int R_sol_size[3];
    float T_sol_data[30];
    int T_sol_size[2];
    //p35p_solver_initialize();


    for(int curr_it = 0; curr_it < it_num; ++curr_it) {
        generateData(X, x, y, f_gen, R_gen, C_gen);

        p35p_single::p35p_solver_single(X, x, y, eps, &solution_num, f_sol_data,
                                        f_sol_size, R_sol_data, R_sol_size, T_sol_data, T_sol_size);

        //allocate for comparison
        float min_diff = INT_MAX; //todo:ok??????????
        float diff_C = INT_MAX;
        float diff_R = INT_MAX;


        for (int i = 0; i < solution_num; ++i) { //todo: make some kind of round
            float f = f_sol_data[i];
            Vector3f T;
            Matrix3f R;
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

            Matrix3f K;
            K.setZero();
            K(0, 0) = K(1, 1) = f;
            K(2, 2) = 1;

            Vector3f C = (-K*R).colPivHouseholderQr().solve(T);

            float diff = abs(f - f_gen);

            float diffR = (R-R_gen).norm()/3;
            float diffC = (C - C_gen).norm()/C_gen.norm();
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
    p35p_single::p35p_solver_single_terminate();
    ASSERT_GT(float(succ_num) / (it_num - zero_solutions_num) , pass_threshold);
}


TEST(AccuracyCheck, P4PDouble){
    double eps = 1e-5;
    int it_num = 10000; //iterations
    int succ_num = 0;
    double succ_threshold = 1.5e-9;
    double pass_threshold = 0.9;
    int zero_solutions_num = 0;

    //allocate for generated data
    double X[12];
    double x[4];
    double y[4];
    double f_gen; Matrix3d R_gen; Vector3d C_gen;

    //allocate for estimated data
    double solution_num;
    double f_sol_data[10];
    int f_sol_size[2];
    double R_sol_data[90];
    int R_sol_size[3];
    double T_sol_data[30];
    int T_sol_size[2];
    //p4p_solver_initialize();


    for(int curr_it = 0; curr_it < it_num; ++curr_it) {
        generateData(X, x, y, f_gen, R_gen, C_gen);

        p4p_double::solve_P4Pf_double(X, x, y, eps, &solution_num, f_sol_data,
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

            Vector3d C = (-R).colPivHouseholderQr().solve(T);

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
        if(min_diff / f_gen < succ_threshold && solution_num != 0)
            succ_num++;
        if(solution_num == 0)
            zero_solutions_num++;
    }
    p4p_double::solve_P4Pf_double_terminate();
    ASSERT_GT(double(succ_num) / (it_num - zero_solutions_num) , pass_threshold);
    //todo: add check for T and C
}


TEST(AccuracyCheck, P4PMixed) {
    float eps = 1e-5;
    int it_num = 10000; //iterations
    int succ_num = 0;
    float succ_threshold = 0.009;
    float pass_threshold = 0.9;
    int zero_solutions_num = 0;

    //allocate for generated data
    float X[12];
    float x[4];
    float y[4];
    float f_gen; Matrix3f R_gen; Vector3f C_gen;

    //allocate for estimated data
    float solution_num;
    float f_sol_data[10];
    int f_sol_size[2];
    float R_sol_data[90];
    int R_sol_size[3];
    float T_sol_data[30];
    int T_sol_size[2];
    //p4p_solver_initialize();


    for(int curr_it = 0; curr_it < it_num; ++curr_it) {
        generateData(X, x, y, f_gen, R_gen, C_gen);

        p4p_mixed::solve_P4Pf_mixed(X, x, y, eps, &solution_num, f_sol_data,
                                    f_sol_size, R_sol_data, R_sol_size, T_sol_data, T_sol_size);

        //allocate for comparison
        float min_diff = INT_MAX; //todo:ok??????????
        float diff_C = INT_MAX;
        float diff_R = INT_MAX;


        for (int i = 0; i < solution_num; ++i) { //todo: make some kind of round
            float f = f_sol_data[i];
            Vector3f T;
            Matrix3f R;
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

            Matrix3f K;
            K.setZero();
            K(0, 0) = K(1, 1) = f;
            K(2, 2) = 1;

            Vector3f C = (-R).colPivHouseholderQr().solve(T);

            float diff = abs(f - f_gen);

            float diffR = (R-R_gen).norm()/3;
            float diffC = (C - C_gen).norm()/C_gen.norm();
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
    p4p_mixed::solve_P4Pf_mixed_terminate();
    ASSERT_GT(float(succ_num) / (it_num - zero_solutions_num) , pass_threshold);
    //todo: add check for T and C
}