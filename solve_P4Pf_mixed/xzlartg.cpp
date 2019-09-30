/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlartg.cpp
 *
 * Code generation for function 'xzlartg'
 *
 */

/* Include files */
#include <cmath>
#include "rt_nonfinite.h"
#include "solve_P4Pf.h"
#include "xzlartg.h"
#include "xgeqp3.h"
#include "solve_P4Pf_rtwutil.h"

/* Function Definitions */
void b_xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T *sn)
{
  float y_tmp;
  float scale;
  float b_y_tmp;
  float f2s;
  float f2;
  float fs_re;
  float fs_im;
  float gs_re;
  float gs_im;
  boolean_T guard1 = false;
  float g2s;
  y_tmp = std::abs(f.re);
  scale = y_tmp;
  b_y_tmp = std::abs(f.im);
  if (b_y_tmp > y_tmp) {
    scale = b_y_tmp;
  }

  f2s = std::abs(g.re);
  f2 = std::abs(g.im);
  if (f2 > f2s) {
    f2s = f2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 5.49755814E+11F) {
    do {
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    guard1 = true;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
    } else {
      do {
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    scale = gs_re * gs_re + gs_im * gs_im;
    f2s = scale;
    if (1.0F > scale) {
      f2s = 1.0F;
    }

    if (f2 <= f2s * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        scale = rt_hypotf_snf(gs_re, gs_im);
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = std::sqrt(scale);
        *cs = rt_hypotf_snf(fs_re, fs_im) / g2s;
        if (b_y_tmp > y_tmp) {
          y_tmp = b_y_tmp;
        }

        if (y_tmp > 1.0F) {
          scale = rt_hypotf_snf(f.re, f.im);
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2s = 5.49755814E+11F * f.re;
          f2 = 5.49755814E+11F * f.im;
          scale = rt_hypotf_snf(f2s, f2);
          fs_re = f2s / scale;
          fs_im = f2 / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = std::sqrt(1.0F + scale / f2);
      *cs = 1.0F / f2s;
      scale += f2;
      fs_re = f2s * fs_re / scale;
      fs_im = f2s * fs_im / scale;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

void xzlartg(const creal32_T f, const creal32_T g, float *cs, creal32_T *sn,
             creal32_T *r)
{
  float y_tmp;
  float scale;
  float b_y_tmp;
  float f2s;
  float f2;
  float fs_re;
  float fs_im;
  float gs_re;
  float gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  float g2s;
  y_tmp = std::abs(f.re);
  scale = y_tmp;
  b_y_tmp = std::abs(f.im);
  if (b_y_tmp > y_tmp) {
    scale = b_y_tmp;
  }

  f2s = std::abs(g.re);
  f2 = std::abs(g.im);
  if (f2 > f2s) {
    f2s = f2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 5.49755814E+11F) {
    do {
      count++;
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    f2 = fs_re * fs_re + fs_im * fs_im;
    scale = gs_re * gs_re + gs_im * gs_im;
    f2s = scale;
    if (1.0F > scale) {
      f2s = 1.0F;
    }

    if (f2 <= f2s * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        r->re = rt_hypotf_snf(g.re, g.im);
        r->im = 0.0F;
        scale = rt_hypotf_snf(gs_re, gs_im);
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = std::sqrt(scale);
        *cs = rt_hypotf_snf(fs_re, fs_im) / g2s;
        if (b_y_tmp > y_tmp) {
          y_tmp = b_y_tmp;
        }

        if (y_tmp > 1.0F) {
          scale = rt_hypotf_snf(f.re, f.im);
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2 = 5.49755814E+11F * f.re;
          f2s = 5.49755814E+11F * f.im;
          scale = rt_hypotf_snf(f2, f2s);
          fs_re = f2 / scale;
          fs_im = f2s / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = std::sqrt(1.0F + scale / f2);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0F / f2s;
      scale += f2;
      f2 = r->re / scale;
      f2s = r->im / scale;
      sn->re = f2 * gs_re - f2s * -gs_im;
      sn->im = f2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 5.49755814E+11F;
          r->im *= 5.49755814E+11F;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 0; rescaledir <= count; rescaledir++) {
            r->re *= 1.8189894E-12F;
            r->im *= 1.8189894E-12F;
          }
        }
      }
    }
  }
}

/* End of code generation (xzlartg.cpp) */
