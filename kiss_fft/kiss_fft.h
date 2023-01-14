#ifndef KISS_FFT_H
#define KISS_FFT_H

#define NFFT 512
#include <math.h>
#include <stddef.h>
#define MAXFACTORS 8
/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned
 4*4*4*2
 */


#define S_MUL(a,b) ( (a)*(b) )
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#define C_MULC(m,a,b) \
    do{ (m).r = (a).r*(b).r + (a).i*(b).i;\
        (m).i = (a).i*(b).r - (a).r*(b).i; }while(0)

#define C_MULBYSCALAR( c, s ) \
    do{ (c).r *= (s);\
        (c).i *= (s); }while(0)

#define  C_ADD( res, a, b) \
    do { (res).r=(a).r+(b).r;\
         (res).i=(a).i+(b).i; }while(0)
#define  C_SUB( res, a, b)\
    do { (res).r=(a).r-(b).r;\
         (res).i=(a).i-(b).i; }while(0)
#define C_ADDTO( res, a)\
    do { (res).r += (a).r;\
         (res).i += (a).i; }while(0)
#define C_SUBFROM( res, a)\
    do { (res).r -= (a).r;\
         (res).i -= (a).i; }while(0)

#  define HALF_OF(x) ((x)*0.5f)

#define  kf_cexp(x,phase) \
        do{ (x)->r = cosf(phase);\
            (x)->i = sinf(phase); }while(0)

typedef struct {
    float r;
    float i;
} cfloat;

typedef struct arch_fft_state{
   int is_supported;
   void *priv;
} arch_fft_state;

typedef struct kiss_fft_state{
    int nfft;
    float scale;
    int shift;
    short factors[2*MAXFACTORS];
    const short *bitrev;
    const cfloat *twiddles;
    arch_fft_state *arch_fft;
} kiss_fft_state;

kiss_fft_state *opus_fft_alloc_twiddles(void * mem, size_t * lenmem, const kiss_fft_state *base, int arch);

kiss_fft_state *opus_fft_alloc(void * mem, size_t * lenmem, int arch);

/**
 * opus_fft(cfg,in_out_buf)
 *
 * Perform an FFT on a complex input buffer.
 * for a forward FFT,
 * fin should be  f[0] , f[1] , ... ,f[nfft-1]
 * fout will be   F[0] , F[1] , ... ,F[nfft-1]
 * Note that each element is complex and can be accessed like
    f[k].r and f[k].i
 * */
void opus_fft_c(const kiss_fft_state *cfg,const cfloat *fin,cfloat *fout);
void opus_ifft_c(const kiss_fft_state *cfg,const cfloat *fin,cfloat *fout);

void opus_fft_impl(const kiss_fft_state *st,cfloat *fout);
int opus_fft_alloc_arch_c(kiss_fft_state *st);



#define opus_fft_alloc_arch(_st, arch) \
         ((void)(arch), opus_fft_alloc_arch_c(_st))

#define opus_fft(_cfg, _fin, _fout, arch) \
         ((void)(arch), opus_fft_c(_cfg, _fin, _fout))

#define opus_ifft(_cfg, _fin, _fout, arch) \
         ((void)(arch), opus_ifft_c(_cfg, _fin, _fout))

#endif
