/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdlanv2.cpp
 *
 * Code generation for function 'xdlanv2'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "p35p_solver.h"
#include "xdlanv2.h"
#include "p35p_solver_rtwutil.h"

/* Function Definitions */
void xdlanv2(float *a, float *b, float *c, float *d, float *rt1r, float *rt1i,
             float *rt2r, float *rt2i, float *cs, float *sn)
{
  float f2;
  float p;
  float bcmis;
  float scale;
  float bcmax;
  float b_scale;
  int b_b;
  int b_c;
  float z;
  float tau;
  float b_p;
  int c_scale;
  if (*c == 0.0F) {
    *cs = 1.0F;
    *sn = 0.0F;
  } else if (*b == 0.0F) {
    *cs = 0.0F;
    *sn = 1.0F;
    bcmis = *d;
    *d = *a;
    *a = bcmis;
    *b = -*c;
    *c = 0.0F;
  } else {
    f2 = *a - *d;
    if ((f2 == 0.0F) && ((*b < 0.0F) != (*c < 0.0F))) {
      *cs = 1.0F;
      *sn = 0.0F;
    } else {
      p = 0.5F * f2;
      scale = std::abs(*b);
      bcmis = std::abs(*c);
      if ((scale > bcmis) || rtIsNaNF(bcmis)) {
        bcmax = scale;
      } else {
        bcmax = bcmis;
      }

      if ((scale < bcmis) || rtIsNaNF(bcmis)) {
        b_scale = scale;
      } else {
        b_scale = bcmis;
      }

      if (!(*b < 0.0F)) {
        b_b = 1;
      } else {
        b_b = -1;
      }

      if (!(*c < 0.0F)) {
        b_c = 1;
      } else {
        b_c = -1;
      }

      bcmis = b_scale * static_cast<float>(b_b) * static_cast<float>(b_c);
      scale = std::abs(p);
      if ((scale > bcmax) || rtIsNaNF(bcmax)) {
      } else {
        scale = bcmax;
      }

      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817842E-16F) {
        *a = std::sqrt(scale) * std::sqrt(z);
        if (!(p < 0.0F)) {
          b_p = *a;
        } else {
          b_p = -*a;
        }

        z = p + b_p;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = rt_hypotf_snf(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0F;
      } else {
        scale = *b + *c;
        tau = rt_hypotf_snf(scale, f2);
        *cs = std::sqrt(0.5F * (1.0F + std::abs(scale) / tau));
        if (!(scale < 0.0F)) {
          c_scale = 1;
        } else {
          c_scale = -1;
        }

        *sn = -(p / (tau * *cs)) * static_cast<float>(c_scale);
        bcmax = *a * *cs + *b * *sn;
        bcmis = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        scale = -*c * *sn + *d * *cs;
        *b = bcmis * *cs + scale * *sn;
        *c = -bcmax * *sn + z * *cs;
        bcmis = 0.5F * ((bcmax * *cs + z * *sn) + (-bcmis * *sn + scale * *cs));
        *a = bcmis;
        *d = bcmis;
        if (*c != 0.0F) {
          if (*b != 0.0F) {
            if ((*b < 0.0F) == (*c < 0.0F)) {
              f2 = std::sqrt(std::abs(*b));
              scale = std::sqrt(std::abs(*c));
              *a = f2 * scale;
              if (!(*c < 0.0F)) {
                p = *a;
              } else {
                p = -*a;
              }

              tau = 1.0F / std::sqrt(std::abs(*b + *c));
              *a = bcmis + p;
              *d = bcmis - p;
              *b -= *c;
              *c = 0.0F;
              bcmax = f2 * tau;
              scale *= tau;
              bcmis = *cs * bcmax - *sn * scale;
              *sn = *cs * scale + *sn * bcmax;
              *cs = bcmis;
            }
          } else {
            *b = -*c;
            *c = 0.0F;
            bcmis = *cs;
            *cs = -*sn;
            *sn = bcmis;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0F) {
    *rt1i = 0.0F;
    *rt2i = 0.0F;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
  }
}

/* End of code generation (xdlanv2.cpp) */
