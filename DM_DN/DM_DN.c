#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "DM_DN.h"

#define USING
CommonState common;
/* Block states (default storage) */
DW_DM_DN_T DM_DN_DW;

/* External inputs (root inport signals with default storage) */
ExtU_DM_DN_T DM_DN_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_DM_DN_T DM_DN_Y;

cfloat N_S[FRAME_SIZE << 1];
cfloat S_N_S[FRAME_SIZE << 1];
cfloat sig[FRAME_SIZE << 1];
cfloat b_x[FRAME_SIZE + 1];
float noisy_s[FRAME_SIZE << 1];
float se_noise_s[FRAME_SIZE << 1];
float sig_m[FRAME_SIZE << 1];
float noisy_ps[FRAME_SIZE + 1];
float se_noise_ps[FRAME_SIZE + 1];
float gammak[FRAME_SIZE + 1];
float ksi[FRAME_SIZE + 1];
float noise_ps[FRAME_SIZE + 1];              /* '<S1>/MATLAB Function1' */
float qk[FRAME_SIZE + 1];                    /* '<S1>/MATLAB Function1' */
float ksi_old[FRAME_SIZE + 1];               /* '<S1>/MATLAB Function1' */
float Xk_prev[FRAME_SIZE + 1];               /* '<S1>/MATLAB Function1' */
float Pglob[FRAME_SIZE + 1];
float b_y1[FRAME_SIZE + 1];
float x2[FRAME_SIZE + 1];
float y2[FRAME_SIZE + 1];
float x2_c[FRAME_SIZE + 1];
float y2_k[FRAME_SIZE + 1];
float save_prev[FRAME_SIZE];             /* '<S1>/MATLAB Function1' */



float rt_hypotf(float u0, float u1);
static void check_init();
static void apply_window(float* x, int N);
static void forward_transform(cfloat* out, const float* in);
static void inverse_transform(float* out, const cfloat* in);
static void DM_DN_expint(const float x[FRAME_SIZE + 1], cfloat y[FRAME_SIZE + 1]);
static void DM_DN_smoothing(const float x[FRAME_SIZE + 1], float y[FRAME_SIZE + 1]);
static void DM_DN_est_sap(float qk[FRAME_SIZE + 1], const float xsi_old[FRAME_SIZE + 1],
                          float zetak[FRAME_SIZE + 1], float *zeta_fr_old, float *z_peak);

float rt_hypotf(float u0, float u1)
{
    float a,b,y;
    a = fabsf(u0);
    b = fabsf(u1);
    if (a < b) {
        a /= b;
        y = sqrtf(a * a + 1.0F) * b;
    }
    else if (a > b) {
        b /= a;
        y = sqrtf(b * b + 1.0F) * a;
    }
    else {
        y = a * 1.41421354F;
    }
    return y;
}

static void check_init() 
{
    int i;
    common.kfft = opus_fft_alloc_twiddles(NULL, NULL, NULL, 0);
    for (i = 0; i < FRAME_SIZE; i++) 
        common.half_window[i] = sinf(.5F * PI * sinf(.5F * PI * (i + .5F) / FRAME_SIZE) * sinf(.5F * PI * (i + .5F) / FRAME_SIZE));
    //for (i = 0; i < FRAME_SIZE; i++)
    //    common.weight[i] = 1.0F/(common.half_window[i] * common.half_window[i] + common.half_window[FRAME_SIZE - 1 - i] * common.half_window[FRAME_SIZE - 1 - i]);

    float norm_c = 0.0F;
    float factor = 1.0F / (float)((SMOOTH_WIN_SIZE << 1) + 1);
    for (i = 0; i < SMOOTH_WIN_SIZE; i++) {
        common.half_smooth_window[i] = sinf(.5F * PI * sinf(.5F * PI * factor * (2 * i + 1)) * sinf(.5F * PI * factor * (2 * i + 1)));
        norm_c += common.half_smooth_window[i] * common.half_smooth_window[i];     
        //printf("win is %f \n", common.half_smooth_window[i]);
    }
    common.smooth_weight = 1.0F / sqrtf(norm_c*2 + 1.0F);
}


static void apply_window(float* x, int N) {
    if (N % 2 == 0)
        for (int i = 0; i < N; i++) {
            x[i] *= common.half_window[i];
            x[(N << 1) - 1 - i] *= common.half_window[i];
        }
    else {
        for (int i = 0; i < N; i++) {
            x[i] *= common.half_window[i];
            x[(N << 1) - i] *= common.half_window[i];
        }
    }
}

static void forward_transform(cfloat* out, const float* in) 
{
    int i;
    cfloat x[FRAME_SIZE << 1];
    for (i = 0; i < FRAME_SIZE << 1; i++) {
        x[i].r = in[i];
        x[i].i = 0.0F;
    }
    opus_fft(common.kfft, x, out, 0);
}

static void inverse_transform(float* out, const cfloat* in) 
{
    int i;
    cfloat x[FRAME_SIZE << 1];
    cfloat y[FRAME_SIZE << 1];
    for (i = 0; i < FRAME_SIZE + 1; i++) {
        x[i] = in[i];
    }
    for (i = FRAME_SIZE + 1; i < (FRAME_SIZE << 1); i++) {
        x[i].r = x[(FRAME_SIZE << 1) - i].r;
        x[i].i = -x[(FRAME_SIZE << 1) - i].i;
    }
    opus_fft(common.kfft, x, y, 0);
    /* output in reverse order for IFFT. */
    out[0] = (FRAME_SIZE << 1) * y[0].r;
    for (i = 1; i < (FRAME_SIZE << 1); i++) 
        out[i] = (FRAME_SIZE << 1) * y[(FRAME_SIZE << 1) - i].r;
}


static void DM_DN_expint(const float x[FRAME_SIZE + 1], cfloat y[FRAME_SIZE + 1])
{
  int k,b_k,exponent, b_exponent,exitg1;
  float j, alpha1, pk, b_x_cv, d, s;
  float am1_im, am1_re, am2_im, am2_re;
  float b_im, bm1_re, bm2_im, bm2_re;
  float f_im, f_re, oldf_im, oldf_re;
  static const float c[9] = { -3.6026937E-9F, -4.81953862E-7F,
    -2.56949825E-5F, -0.000697379059F, -0.0101957349F, -0.0781186372F,
    -0.301243275F, -0.777380705F, 8.26766205F };

  for (k = 0; k < FRAME_SIZE + 1; k++) {
    b_x_cv = x[k];
    pk = -3.6026937E-9F;
    for (b_k = 0; b_k < 8; b_k++) 
      pk = b_x_cv * pk + c[b_k + 1];

    if (x[k] == 0.0F) {
      oldf_re = 3.402823466E+38F;
      oldf_im = 0.0F;
    } 
    else if (pk >= 0.0F) {
      oldf_re = x[k];
      oldf_im = 0.0F;
      if (oldf_re < 0.0F) {
        oldf_re = logf(fabsf(oldf_re));
        oldf_im = PI;
      } 
      else 
        oldf_re = logf(oldf_re);


      oldf_re = -0.577215672F - oldf_re;
      oldf_im = 0.0F - oldf_im;
      j = 1.0F;
      pk = x[k];
      f_re = x[k];
      do {
        exitg1 = 0;
        f_im = rt_hypotf(oldf_re, oldf_im);
        if (f_im <= 1.17549435E-38F) {
          f_im = 1.4013E-45F;
        } else {
          frexpf(f_im, &exponent);
          f_im = ldexpf(1.0F, exponent - 24);
        }

        if (fabsf(f_re) > f_im) {
          oldf_re += f_re;
          j += 1.0F;
          pk = -b_x_cv * pk / j;
          f_re = pk / j;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    } 
    else if (pk < 0.0F) {
      am2_re = 0.0F;
      am2_im = 0.0F;
      bm2_re = 1.0F;
      bm2_im = 0.0F;
      am1_re = 1.0F;
      am1_im = 0.0F;
      bm1_re = x[k];
      f_re = 1.0F / bm1_re;
      f_im = 0.0F;
      oldf_re = 3.402823466E+38F;
      oldf_im = 0.0F;
      j = 2.0F;
      do {
        exitg1 = 0;
        pk = rt_hypotf(f_re, f_im);
        if (pk > 0.0F) {
          frexpf(pk, &b_exponent);
          pk = ldexpf(1.0F, b_exponent - 24);
        }

        oldf_re = f_re - oldf_re;
        oldf_im = f_im - oldf_im;
        if (rt_hypotf(oldf_re, oldf_im) > 100.0F * pk) {
          alpha1 = j * 0.5F;
          oldf_re = alpha1 * am2_re + am1_re;
          oldf_im = alpha1 * am2_im + am1_im;
          pk = alpha1 * bm2_re + bm1_re;
          b_im = alpha1 * bm2_im;
          if (b_im == 0.0F) {
            if (am1_im == 0.0F) {
              am2_re = am1_re / pk;
              am2_im = 0.0F;
            } 
            else if (am1_re == 0.0F) {
              am2_re = 0.0F;
              am2_im = am1_im / pk;
            } 
            else {
              am2_re = am1_re / pk;
              am2_im = am1_im / pk;
            }
          } 
          else if (pk == 0.0F) {
            if (am1_re == 0.0F) {
              am2_re = am1_im / b_im;
              am2_im = 0.0F;
            } 
            else if (am1_im == 0.0F) {
              am2_re = 0.0F;
              am2_im = -(am1_re / b_im);
            } 
            else {
              am2_re = am1_im / b_im;
              am2_im = -(am1_re / b_im);
            }
          } 
          else {
            s = pk / b_im;
            d = s * pk + b_im;
            am2_re = (s * am1_re + am1_im) / d;
            am2_im = (s * am1_im - am1_re) / d;
          }

          if (b_im == 0.0F) {
            if (oldf_im == 0.0F) {
              am1_re = oldf_re / pk;
              am1_im = 0.0F;
            } 
            else if (oldf_re == 0.0F) {
              am1_re = 0.0F;
              am1_im = oldf_im / pk;
            } 
            else {
              am1_re = oldf_re / pk;
              am1_im = oldf_im / pk;
            }
          } 
          else if (pk == 0.0F) {
            if (oldf_re == 0.0F) {
              am1_re = oldf_im / b_im;
              am1_im = 0.0F;
            } 
            else if (oldf_im == 0.0F) {
              am1_re = 0.0F;
              am1_im = -(oldf_re / b_im);
            } 
            else {
              am1_re = oldf_im / b_im;
              am1_im = -(oldf_re / b_im);
            }
          } 
          else {
            s = pk / b_im;
            d = s * pk + b_im;
            am1_re = (s * oldf_re + oldf_im) / d;
            am1_im = (s * oldf_im - oldf_re) / d;
          }

          f_re = am1_re;
          f_im = am1_im;
          alpha1 = ((j + 1.0F) - 1.0F) / 2.0F;
          oldf_re = b_x_cv * am1_re + alpha1 * am2_re;
          oldf_im = b_x_cv * am1_im + alpha1 * am2_im;
          if (b_im == 0.0F) {
            am2_re = bm1_re / pk;
            bm1_re = 0.0F;
          } 
          else if (pk == 0.0F) {
            if (bm1_re == 0.0F) {
              am2_re = 0.0F / b_im;
              bm1_re = 0.0F;
            } 
            else {
              am2_re = 0.0F;
              bm1_re = -(bm1_re / b_im);
            }
          } 
          else {
            s = pk / b_im;
            d = s * pk + b_im;
            am2_re = s * bm1_re / d;
            bm1_re = (0.0F - bm1_re) / d;
          }

          pk = alpha1 * am2_re + b_x_cv;
          b_im = alpha1 * bm1_re;
          if (b_im == 0.0F) {
            if (am1_im == 0.0F) {
              am2_re = am1_re / pk;
              am2_im = 0.0F;
            } 
            else if (am1_re == 0.0F) {
              am2_re = 0.0F;
              am2_im = am1_im / pk;
            } 
            else {
              am2_re = am1_re / pk;
              am2_im = am1_im / pk;
            }
          } 
          else if (pk == 0.0F) {
            if (am1_re == 0.0F) {
              am2_re = am1_im / b_im;
              am2_im = 0.0F;
            } 
            else if (am1_im == 0.0F) {
              am2_re = 0.0F;
              am2_im = -(am1_re / b_im);
            } 
            else {
              am2_re = am1_im / b_im;
              am2_im = -(am1_re / b_im);
            }
          } 
          else {
            s = pk / b_im;
            d = s * pk + b_im;
            am2_re = (s * am1_re + am1_im) / d;
            am2_im = (s * am1_im - am1_re) / d;
          }

          if (b_im == 0.0F) {
            bm2_re = 1.0F / pk;
            bm2_im = 0.0F;
          } 
          else if (pk == 0.0F) {
            bm2_re = 0.0F;
            bm2_im = -(1.0F / b_im);
          } 
          else {
            s = pk / b_im;
            d = s * pk + b_im;
            bm2_re = s / d;
            bm2_im = -1.0F / d;
          }

          if (b_im == 0.0F) {
            if (oldf_im == 0.0F) {
              am1_re = oldf_re / pk;
              am1_im = 0.0F;
            } 
            else if (oldf_re == 0.0F) {
              am1_re = 0.0F;
              am1_im = oldf_im / pk;
            } 
            else {
              am1_re = oldf_re / pk;
              am1_im = oldf_im / pk;
            }
          } 
          else if (pk == 0.0F) {
            if (oldf_re == 0.0F) {
              am1_re = oldf_im / b_im;
              am1_im = 0.0F;
            } 
            else if (oldf_im == 0.0F) {
              am1_re = 0.0F;
              am1_im = -(oldf_re / b_im);
            } 
            else {
              am1_re = oldf_im / b_im;
              am1_im = -(oldf_re / b_im);
            }
          } 
          else {
            s = pk / b_im;
            d = s * pk + b_im;
            am1_re = (s * oldf_re + oldf_im) / d;
            am1_im = (s * oldf_im - oldf_re) / d;
          }

          bm1_re = 1.0F;
          oldf_re = f_re;
          oldf_im = f_im;
          f_re = am1_re;
          f_im = am1_im;
          j = (j + 1.0F) + 1.0F;
        } 
        else 
          exitg1 = 1;

      } while (exitg1 == 0);

      b_x_cv = expf(-x[k]);
      oldf_re = b_x_cv * f_re;
      oldf_im = b_x_cv * f_im;
      if (x[k] < 0.0F) 
        oldf_im -= 3.14159274F;
    } else {
      oldf_re = 0.0F;
      oldf_im = 0.0F;
    }

    y[k].r = oldf_re;
    y[k].i = oldf_im;
  }
}

/* Function for MATLAB Function: '<S1>/MATLAB Function1' */
static void DM_DN_smoothing(const float x[FRAME_SIZE + 1], float y[FRAME_SIZE + 1])
{
  int i,j;
  memset(&y[0], 0, (FRAME_SIZE + 1) * sizeof(float));
  for (i = 0; i < SMOOTH_WIN_SIZE + 1; i++)
      for (j = i; j < FRAME_SIZE + 1; j++) {
          if (i != SMOOTH_WIN_SIZE)
              y[j] += x[j - i] * common.half_smooth_window[i];
          else
              y[j] += x[j - i];// common.half_smooth_window[SMOOTH_WIN_SIZE] = 1
      }
  memset(&x2_c[0], 0, (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&x2_c[0], &x[SMOOTH_WIN_SIZE], (FRAME_SIZE + 1 - SMOOTH_WIN_SIZE) * sizeof(float));
  memset(&y2_k[0], 0, (FRAME_SIZE + 1) * sizeof(float));

  for (i = 0; i < SMOOTH_WIN_SIZE; i++)
      for (j = i; j < FRAME_SIZE + 1; j++)
          y2_k[j] += x2_c[j - i] * common.half_smooth_window[SMOOTH_WIN_SIZE - 1 - i];
  
  for (i = 0; i < FRAME_SIZE + 1; i++)
      y[i] = (y[i] + y2_k[i]) * common.smooth_weight;

}

/* Function for MATLAB Function: '<S1>/MATLAB Function1' */
static void DM_DN_est_sap(float qk[FRAME_SIZE + 1], const float xsi_old[FRAME_SIZE + 1],
                          float zetak[FRAME_SIZE + 1], float *zeta_fr_old, float *z_peak)
{
  int k;
  float Pframe, y, zeta_fr;
  for (k = 0; k < FRAME_SIZE + 1; k++) {
    zetak[k] = 0.7F * zetak[k] + 0.3F * xsi_old[k];
    b_y1[k] = 0.382683F * zetak[k];;
  }

  for (k = 1; k < FRAME_SIZE + 1; k++) {
    b_y1[k] += zetak[k - 1];
  }

  memset(&x2[0], 0, (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&x2[0], &zetak[1], FRAME_SIZE * sizeof(float));
  for (k = 0; k < FRAME_SIZE + 1; k++) {
    Pglob[k] = 0.0F;
    y2[k] = 0.0F;
    b_y1[k] = (0.382683F * x2[k] + b_y1[k]) / sqrtf(1.292893F);
  }

  DM_DN_smoothing(zetak, x2);
  for (k = 0; k < FRAME_SIZE + 1; k++) {
    Pframe = x2[k];
    zeta_fr = b_y1[k];
    if (zeta_fr > 0.3162F) 
      y2[k] = 1.0F;

    if ((zeta_fr > 0.1F) && (zeta_fr < 0.3162F)) 
      y2[k] = log10f(zeta_fr * 10.0F) / 0.499961853F;

    if (Pframe > 0.3162F) 
      Pglob[k] = 1.0F;

    if ((Pframe > 0.1F) && (Pframe < 0.3162F)) 
      Pglob[k] = log10f(Pframe * 10.0F) / 0.499961853F;
  }

  zeta_fr = zetak[0];
  for (k = 0; k < FRAME_SIZE; k++) 
    zeta_fr += zetak[k + 1];

  zeta_fr /= FRAME_SIZE + 1.0F;
  if (zeta_fr > 0.1F) {
    if (zeta_fr > *zeta_fr_old) {
      Pframe = 1.0F;
      y = fmaxf(zeta_fr, 1.0F);
      y = fminf(y, 10.0F);
    } 
    else if (zeta_fr <= *z_peak * 0.1F)
      Pframe = 0.0F;
    else if (zeta_fr >= *z_peak * 0.3162F)
      Pframe = 1.0F;
    else
      Pframe = log10f(zeta_fr / *z_peak * 10.0F) / 0.499961853F;
  } 
  else
    Pframe = 0.0F;

  *zeta_fr_old = zeta_fr;
  for (k = 0; k < FRAME_SIZE + 1; k++) {
    zeta_fr = 1.0F - y2[k] * Pglob[k] * Pframe;
    qk[k] = fmaxf(zeta_fr, 0.95F);
    y2[k] = zeta_fr;
  }
}

void DM_DN_step(void)
{
  int i;
  float zeta_fr_old, counter, counter_b, ksi_p, r;

  memcpy(&noise_ps[0], &DM_DN_DW.Delay7_DSTATE[0], (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&qk[0], &DM_DN_DW.Delay8_DSTATE[0], (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&ksi_old[0], &DM_DN_DW.Delay9_DSTATE[0], (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&Xk_prev[0], &DM_DN_DW.Delay10_DSTATE[0], (FRAME_SIZE + 1) * sizeof(float));

  zeta_fr_old = DM_DN_DW.Delay12_DSTATE;


  for (i = 0; i < FRAME_SIZE; i++) {
    save_prev[i] = DM_DN_DW.Delay14_DSTATE[i];
    noisy_s[i] = DM_DN_DW.Delay_DSTATE[i];
    noisy_s[i + FRAME_SIZE] = DM_DN_U.noisy_speech[i];
  }

    
  counter = DM_DN_DW.Delay1_DSTATE + 1.0F;
  if (DM_DN_DW.Delay1_DSTATE + 1.0F > DM_DN_U.NIS) 
    counter = DM_DN_U.NIS + 1.0F;

  apply_window(noisy_s, FRAME_SIZE);
  forward_transform(N_S, noisy_s);

  for (i = 0; i < FRAME_SIZE + 1; i++) {
    counter_b = rt_hypotf(N_S[i].r,N_S[i].i);
    noisy_ps[i] = counter_b * counter_b;
  }

  if (counter < DM_DN_U.NIS) {
    memcpy(&DM_DN_Y.clean_speech[0], &noisy_s[0], FRAME_SIZE * sizeof(float));
    for (i = 0; i < FRAME_SIZE + 1; i++) {
      counter_b = noise_ps[i] + noisy_ps[i];
      if (DM_DN_U.NIS - 1.0F == counter) 
        counter_b /= DM_DN_U.NIS - 1.0F;

      noise_ps[i] = counter_b;
    }
  } 
  else {
    for (i = 0; i < FRAME_SIZE + 1; i++) {
      counter_b = noisy_ps[i] / noise_ps[i];
      gammak[i] = fminf(counter_b, 40.0F);
    }

    if (counter == DM_DN_U.NIS) {
      for (i = 0; i < FRAME_SIZE + 1; i++) {
        counter_b = gammak[i] - 1.0F;
        if (counter_b <= 0.0F) 
          counter_b = 0.0F;

        ksi[i] = 0.02F * counter_b + 0.98F;
        qk[i] += 0.5F;
      }

      zeta_fr_old = DM_DN_DW.Delay12_DSTATE + 1000.0F;
    } 
    else {
      for (i = 0; i < FRAME_SIZE + 1; i++) {
        counter_b = gammak[i] - 1.0F;
        if (counter_b <= 0.0F) 
          counter_b = 0.0F;

        counter_b = 0.98F * Xk_prev[i] / noise_ps[i] + 0.02F * counter_b;
        ksi[i] = fmaxf(counter_b, 0.00316227763F);
      }
    }

    for (i = 0; i < FRAME_SIZE + 1; i++) {
      ksi_p = ksi[i];
      ksi_p /= ksi_p + 1.0F;
      Xk_prev[i] = ksi_p;
      gammak[i] *= ksi_p;
    }

    DM_DN_expint(gammak, b_x);
    for (i = 0; i < FRAME_SIZE + 1; i++) {
      if (b_x[i].i == 0.0F) {
        ksi_p = expf(0.5F * b_x[i].r);
        r = 0.0F;
      } 
      else {
        r = expf(0.5F * b_x[i].r / 2.0F);
        ksi_p = r * cosf(0.5F * b_x[i].i) * r;
        r *= r * sinf(0.5F * b_x[i].i);
      }

      b_x[i].r = ksi_p;
      b_x[i].i = r;
      Xk_prev[i] *= ksi_p;
    }

    if (DM_DN_U.gate) {
      DM_DN_est_sap(qk, ksi_old, DM_DN_DW.Delay11_DSTATE,&zeta_fr_old, &DM_DN_DW.Delay13_DSTATE);
      for (i = 0; i < FRAME_SIZE + 1; i++) {
        ksi_p = qk[i];
        ksi_p = (1.0F - ksi_p) / ((ksi[i] + 1.0F) * ksi_p * expf(-gammak[i]) + (1.0F - ksi_p));         
        gammak[i] = ksi_p;
        Xk_prev[i] = powf(0.01F, 1.0F - ksi_p) * powf(Xk_prev[i], ksi_p);
      }
    }

    for (i = 0; i < FRAME_SIZE << 1; i++) 
      noisy_s[i] = rt_hypotf(N_S[i].r, N_S[i].i);
    
    for (i = 0; i < FRAME_SIZE + 1; i++) 
      sig_m[i] = noisy_s[i] * Xk_prev[i];

    for (i = 0; i < FRAME_SIZE - 1; i++)
      sig_m[i + FRAME_SIZE + 1] = noisy_s[i + FRAME_SIZE + 1] * Xk_prev[FRAME_SIZE - 1 - i];

    for (i = 0; i < FRAME_SIZE + 1; i++) {
      Xk_prev[i] = sig_m[i] * sig_m[i];
      ksi_old[i] = ksi[i];
    }

    for (i = 0; i < FRAME_SIZE << 1; i++) {
      r = atan2f(N_S[i].i, N_S[i].r);
      ksi_p = cosf(r);
      r = sinf(r);

      sig[i].r = sig_m[i] * ksi_p;
      sig[i].i = sig_m[i] * r;
      N_S[i].r = ksi_p;
      N_S[i].i = r;
    }
   
    inverse_transform(sig_m, sig);
    apply_window(sig_m, FRAME_SIZE);
    for (i = 0; i < (FRAME_SIZE << 1); i++)
        sig_m[i] = fmaxf(fminf(sig_m[i], 1.0F), -1.0F);
    for (i = 0; i < FRAME_SIZE; i++) {
        DM_DN_Y.clean_speech[i] = (sig_m[i] + save_prev[i]);// *common.weight[i];  // for i in common.weight == 1 in this window   
        save_prev[i] = sig_m[i + FRAME_SIZE];
    }

    for (i = 0; i < FRAME_SIZE; i++) {
        se_noise_s[i] = DM_DN_DW.Delay2_DSTATE[i];
        se_noise_s[i + FRAME_SIZE] = DM_DN_U.se_noise[i];
    }

    apply_window(se_noise_s, FRAME_SIZE);
    forward_transform(S_N_S, se_noise_s);

    for (i = 0; i < FRAME_SIZE + 1; i++) {
        counter_b = rt_hypotf(S_N_S[i].r, S_N_S[i].i);
        se_noise_ps[i] = counter_b * counter_b;
        noise_ps[i] = DM_DN_U.mu * noise_ps[i] + (1.0F - DM_DN_U.mu) * se_noise_ps[i];
    }


            
  }

  memcpy(&DM_DN_DW.Delay_DSTATE[0], &DM_DN_U.noisy_speech[0], FRAME_SIZE * sizeof(float));
  DM_DN_DW.Delay1_DSTATE = counter;
  memcpy(&DM_DN_DW.Delay2_DSTATE[0], &DM_DN_U.se_noise[0], FRAME_SIZE * sizeof(float));
  memcpy(&DM_DN_DW.Delay7_DSTATE[0], &noise_ps[0], (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&DM_DN_DW.Delay8_DSTATE[0], &qk[0], (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&DM_DN_DW.Delay9_DSTATE[0], &ksi_old[0], (FRAME_SIZE + 1) * sizeof(float));
  memcpy(&DM_DN_DW.Delay10_DSTATE[0], &Xk_prev[0], (FRAME_SIZE + 1) * sizeof(float));
  DM_DN_DW.Delay12_DSTATE = zeta_fr_old;
  memcpy(&DM_DN_DW.Delay14_DSTATE[0], &save_prev[0], FRAME_SIZE * sizeof(float));
}


void algorithm_run(uint32_t *input, uint32_t *output, uint32_t len)
{
    uint8_t j;
    
    for (j = 0; j < FRAME_SIZE; j++) {
        
        DM_DN_U.noisy_speech[j] = ((int32_t)input[2*j] << 8) / 2147483648.0f;
    }
    DM_DN_step();
    for (j = 0; j < FRAME_SIZE; j++) {
        output[j*2]   =   (int32_t)(DM_DN_Y.clean_speech[j] * 2147483648.0f) >> 8;;
        output[j*2+1] = output[j*2];
    }

	//for (j = 0; j < FRAME_SIZE; j++) {
 //      output[j*2]   =   input[2*j];
 //      output[j*2+1] = output[j*2];
 //   }
		
}

uint16_t noice_reduction_set(uint8_t *cfg, uint8_t len)
{  
    return 0;
}

uint16_t noice_reduction_get(uint8_t *cfg)
{
    return 0;
}

/* Model initialize function */
void noice_reduction_init(uint32_t samp_rate)
{
    DM_DN_U.NIS = 20.0F;
    DM_DN_U.mu = 0.8F;
    DM_DN_U.gate = 1U;
    check_init();
}

#ifdef USING
#define MAXDATALEN (16000*10)
#include <time.h>
int main() 
{
    uint32_t rate = 8000;
    noice_reduction_init(rate);

    //读文件 
    float* input = (float*)malloc(MAXDATALEN * sizeof(float));
    FILE* fp_input;
    fp_input = fopen("C:\\Users\\ll\\Desktop\\C DM_RT_DN\\ch01.txt", "r"); //读的音频数据的txt,
    if (fp_input == NULL) {
        printf("file is error.");
        return -1;
    }

    for (int i = 0; i < MAXDATALEN; i++) {
        fscanf(fp_input, "%f", &input[i]);
        //printf("%f\n", input[i]);
    }

    float* noise = (float*)malloc(MAXDATALEN * sizeof(float));
    FILE* fp_noise;
    fp_noise = fopen("C:\\Users\\ll\\Desktop\\C DM_RT_DN\\ch02.txt", "r"); //读的音频数据的txt,
    if (fp_noise == NULL) {
        printf("file is error.");
        return -1;
    }

    for (int i = 0; i < MAXDATALEN; i++) {
        fscanf(fp_noise, "%f", &noise[i]);
        //printf("%f\n", noise[i]);
    }

    FILE* fp_output;
    fp_output = fopen("C:\\Users\\ll\\Desktop\\C DM_RT_DN\\output_ch.txt", "w");

    clock_t start, finish;
    double Total_time;
    start = clock();

    for (long i = 0; i < MAXDATALEN - FRAME_SIZE; i += FRAME_SIZE) {
        for (int j = 0; j < FRAME_SIZE; j++) {
            DM_DN_U.noisy_speech[j] = input[i + j];
            DM_DN_U.se_noise[j] = noise[i + j];
        }
        DM_DN_step();
        for (int j = 0; j < FRAME_SIZE; j++) {
            fprintf(fp_output, "%f\n", DM_DN_Y.clean_speech[j]);
            printf("%f\n", DM_DN_Y.clean_speech[j]);
        }
    }
    free(input);
    fclose(fp_input);
    free(noise);
    fclose(fp_noise);
    fclose(fp_output);

    finish = clock();
    Total_time = (float)(finish - start) / CLOCKS_PER_SEC; //单位换算成秒
    printf("%f seconds\n", Total_time);
    return 0;
}
#endif
