#ifndef DM_DN_H
#define DM_DN_H

#include <math.h>
#include <string.h>
#include "kiss_fft.h"
#include <stdint.h>

#define FRAME_SIZE (NFFT >> 1)
#define PI 3.1415927F
#define SMOOTH_WIN_SIZE 15

/*for XK speech files*/

#define SQUARE(x) ((x)*(x))
#define SMOOTH_BANDS 1

typedef struct {
    kiss_fft_state* kfft;
    float half_window[FRAME_SIZE];
    //float weight[FRAME_SIZE];
    float half_smooth_window[SMOOTH_WIN_SIZE];
    float smooth_weight;
} CommonState;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  float Delay_DSTATE[FRAME_SIZE];          /* '<S1>/Delay' */
  float Delay1_DSTATE;              /* '<S1>/Delay1' */
  float Delay2_DSTATE[FRAME_SIZE];          /* '<S1>/Delay2' */
  float Delay7_DSTATE[FRAME_SIZE + 1];         /* '<S1>/Delay7' */
  float Delay8_DSTATE[FRAME_SIZE + 1];         /* '<S1>/Delay8' */
  float Delay9_DSTATE[FRAME_SIZE + 1];         /* '<S1>/Delay9' */
  float Delay10_DSTATE[FRAME_SIZE + 1];        /* '<S1>/Delay10' */
  float Delay11_DSTATE[FRAME_SIZE + 1];        /* '<S1>/Delay11' */
  float Delay12_DSTATE;             /* '<S1>/Delay12' */
  float Delay13_DSTATE;             /* '<S1>/Delay13' */
  float Delay14_DSTATE[FRAME_SIZE];        /* '<S1>/Delay14' */

} DW_DM_DN_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  float noisy_speech[FRAME_SIZE];          /* '<Root>/noisy_speech' */
  float se_noise[FRAME_SIZE];
  float NIS;                        /* '<Root>/NIS' */
  float mu;                        /* '<Root>/eta' */
  unsigned int gate;                      /* '<Root>/gate' */
} ExtU_DM_DN_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  float clean_speech[FRAME_SIZE];          /* '<Root>/clean_speech' */
} ExtY_DM_DN_T;



extern CommonState common;
/* Block states (default storage) */
extern DW_DM_DN_T DM_DN_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_DM_DN_T DM_DN_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_DM_DN_T DM_DN_Y;

extern void algorithm_run(uint32_t *input, uint32_t *output, uint32_t len);
extern uint16_t noice_reduction_set(uint8_t *cfg, uint8_t len);
extern uint16_t noice_reduction_get(uint8_t *cfg);
extern void noice_reduction_init(uint32_t samp_rate);

/* Model entry point functions */
extern void DM_DN_initialize(void);
extern void DM_DN_step(void);

#endif
