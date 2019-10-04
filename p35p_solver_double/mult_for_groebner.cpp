//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: mult_for_groebner.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 04-Oct-2019 01:44:03
//

// Include Files
#include "mult_for_groebner.h"
#include "p35p_solver.h"
#include "rt_nonfinite.h"
#include <cstring>

// Function Definitions

//
// Arguments    : const double G4[112]
//                double G20[900]
// Return Type  : void
//
void mult_for_groebner(const double G4[112], double G20[900])
{
  int i;
  double G20_tmp;
  double b_G20_tmp;
  double c_G20_tmp;
  double d_G20_tmp;
  double e_G20_tmp;
  double f_G20_tmp;
  double g_G20_tmp;
  double h_G20_tmp;
  double i_G20_tmp;
  double j_G20_tmp;
  double k_G20_tmp;
  double l_G20_tmp;
  double m_G20_tmp;
  double n_G20_tmp;
  double o_G20_tmp;
  double p_G20_tmp;
  double q_G20_tmp;
  double r_G20_tmp;
  double s_G20_tmp;
  double t_G20_tmp;
  double u_G20_tmp;
  double v_G20_tmp;
  double w_G20_tmp;
  double x_G20_tmp;
  double y_G20_tmp;
  double ab_G20_tmp;
  double bb_G20_tmp;
  int b_i;
  std::memset(&G20[0], 0, 900U * sizeof(double));

  // multiply by qx^2
  // multiply by qxqy
  // multiply by qx
  // multiply by qy
  // multiply by 1
  for (i = 0; i < 4; i++) {
    G20[i] = G4[i];
    G20_tmp = G4[i + 4];
    G20[i + 20] = G20_tmp;
    b_G20_tmp = G4[i + 8];
    G20[i + 40] = b_G20_tmp;
    c_G20_tmp = G4[i + 12];
    G20[i + 60] = c_G20_tmp;
    d_G20_tmp = G4[i + 16];
    G20[i + 80] = d_G20_tmp;
    e_G20_tmp = G4[i + 20];
    G20[i + 100] = e_G20_tmp;
    f_G20_tmp = G4[i + 24];
    G20[i + 120] = f_G20_tmp;
    g_G20_tmp = G4[i + 28];
    G20[i + 180] = g_G20_tmp;
    h_G20_tmp = G4[i + 32];
    G20[i + 200] = h_G20_tmp;
    i_G20_tmp = G4[i + 36];
    G20[i + 220] = i_G20_tmp;
    j_G20_tmp = G4[i + 40];
    G20[i + 240] = j_G20_tmp;
    k_G20_tmp = G4[i + 44];
    G20[i + 260] = k_G20_tmp;
    l_G20_tmp = G4[i + 48];
    G20[i + 280] = l_G20_tmp;
    m_G20_tmp = G4[i + 52];
    G20[i + 340] = m_G20_tmp;
    n_G20_tmp = G4[i + 56];
    G20[i + 360] = n_G20_tmp;
    o_G20_tmp = G4[i + 60];
    G20[i + 380] = o_G20_tmp;
    p_G20_tmp = G4[i + 64];
    G20[i + 400] = p_G20_tmp;
    q_G20_tmp = G4[i + 68];
    G20[i + 420] = q_G20_tmp;
    r_G20_tmp = G4[i + 72];
    G20[i + 480] = r_G20_tmp;
    s_G20_tmp = G4[i + 76];
    G20[i + 500] = s_G20_tmp;
    t_G20_tmp = G4[i + 80];
    G20[i + 520] = t_G20_tmp;
    u_G20_tmp = G4[i + 84];
    G20[i + 540] = u_G20_tmp;
    v_G20_tmp = G4[i + 88];
    G20[i + 600] = v_G20_tmp;
    w_G20_tmp = G4[i + 92];
    G20[i + 620] = w_G20_tmp;
    x_G20_tmp = G4[i + 96];
    G20[i + 640] = x_G20_tmp;
    y_G20_tmp = G4[i + 100];
    G20[i + 700] = y_G20_tmp;
    ab_G20_tmp = G4[i + 104];
    G20[i + 720] = ab_G20_tmp;
    bb_G20_tmp = G4[i + 108];
    G20[i + 780] = bb_G20_tmp;
    G20[i + 24] = G4[i];
    G20[i + 44] = G20_tmp;
    G20[i + 64] = b_G20_tmp;
    G20[i + 84] = c_G20_tmp;
    G20[i + 104] = d_G20_tmp;
    G20[i + 124] = e_G20_tmp;
    G20[i + 144] = f_G20_tmp;
    G20[i + 204] = g_G20_tmp;
    G20[i + 224] = h_G20_tmp;
    G20[i + 244] = i_G20_tmp;
    G20[i + 264] = j_G20_tmp;
    G20[i + 284] = k_G20_tmp;
    G20[i + 304] = l_G20_tmp;
    G20[i + 364] = m_G20_tmp;
    G20[i + 384] = n_G20_tmp;
    G20[i + 404] = o_G20_tmp;
    G20[i + 424] = p_G20_tmp;
    G20[i + 444] = q_G20_tmp;
    G20[i + 504] = r_G20_tmp;
    G20[i + 524] = s_G20_tmp;
    G20[i + 544] = t_G20_tmp;
    G20[i + 564] = u_G20_tmp;
    G20[i + 624] = v_G20_tmp;
    G20[i + 644] = w_G20_tmp;
    G20[i + 664] = x_G20_tmp;
    G20[i + 724] = y_G20_tmp;
    G20[i + 744] = ab_G20_tmp;
    G20[i + 804] = bb_G20_tmp;
    G20[i + 188] = G4[i];
    G20[i + 208] = G20_tmp;
    G20[i + 228] = b_G20_tmp;
    G20[i + 248] = c_G20_tmp;
    G20[i + 268] = d_G20_tmp;
    G20[i + 288] = e_G20_tmp;
    G20[i + 308] = f_G20_tmp;
    G20[i + 348] = g_G20_tmp;
    G20[i + 368] = h_G20_tmp;
    G20[i + 388] = i_G20_tmp;
    G20[i + 408] = j_G20_tmp;
    G20[i + 428] = k_G20_tmp;
    G20[i + 448] = l_G20_tmp;
    G20[i + 488] = m_G20_tmp;
    G20[i + 508] = n_G20_tmp;
    G20[i + 528] = o_G20_tmp;
    G20[i + 548] = p_G20_tmp;
    G20[i + 568] = q_G20_tmp;
    G20[i + 608] = r_G20_tmp;
    G20[i + 628] = s_G20_tmp;
    G20[i + 648] = t_G20_tmp;
    G20[i + 668] = u_G20_tmp;
    G20[i + 708] = v_G20_tmp;
    G20[i + 728] = w_G20_tmp;
    G20[i + 748] = x_G20_tmp;
    G20[i + 788] = y_G20_tmp;
    G20[i + 808] = ab_G20_tmp;
    G20[i + 848] = bb_G20_tmp;
    G20[i + 212] = G4[i];
    G20[i + 232] = G20_tmp;
    G20[i + 252] = b_G20_tmp;
    G20[i + 272] = c_G20_tmp;
    G20[i + 292] = d_G20_tmp;
    G20[i + 312] = e_G20_tmp;
    G20[i + 332] = f_G20_tmp;
    G20[i + 372] = g_G20_tmp;
    G20[i + 392] = h_G20_tmp;
    G20[i + 412] = i_G20_tmp;
    G20[i + 432] = j_G20_tmp;
    G20[i + 452] = k_G20_tmp;
    G20[i + 472] = l_G20_tmp;
    G20[i + 512] = m_G20_tmp;
    G20[i + 532] = n_G20_tmp;
    G20[i + 552] = o_G20_tmp;
    G20[i + 572] = p_G20_tmp;
    G20[i + 592] = q_G20_tmp;
    G20[i + 632] = r_G20_tmp;
    G20[i + 652] = s_G20_tmp;
    G20[i + 672] = t_G20_tmp;
    G20[i + 692] = u_G20_tmp;
    G20[i + 732] = v_G20_tmp;
    G20[i + 752] = w_G20_tmp;
    G20[i + 772] = x_G20_tmp;
    G20[i + 812] = y_G20_tmp;
    G20[i + 832] = ab_G20_tmp;
    G20[i + 872] = bb_G20_tmp;
    for (b_i = 0; b_i < 28; b_i++) {
      G20[(i + 20 * (b_i + 17)) + 16] = G4[i + (b_i << 2)];
    }
  }
}

//
// File trailer for mult_for_groebner.cpp
//
// [EOF]
//
