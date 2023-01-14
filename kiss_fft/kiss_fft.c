#include <math.h>
#include <stdio.h>
#include "kiss_fft.h"
#define CUSTOM_MODES
static void kf_bfly2(cfloat* Fout, int m, int N)
{
   cfloat * Fout2;
   int i;
#ifdef CUSTOM_MODES
   if (m==1){
       cfloat t;
       for (i=0;i<N;i++) {        
         Fout2 = Fout + 1;
         t = *Fout2;
         C_SUB( *Fout2 ,  *Fout , t );
         C_ADDTO( *Fout ,  t );
         Fout += 2;
      }
   } else
#endif
   {
      float tw;
      tw = 0.7071067812f;
      /* We know that m==4 here because the radix-2 is just after a radix-4 */
      for (i=0;i<N;i++) {
         cfloat t;
         Fout2 = Fout + 4;
         t = Fout2[0];
         C_SUB( Fout2[0] ,  Fout[0] , t );
         C_ADDTO( Fout[0] ,  t );

         t.r = S_MUL(Fout2[1].r + Fout2[1].i, tw);
         t.i = S_MUL(Fout2[1].i - Fout2[1].r, tw);
         C_SUB( Fout2[1] ,  Fout[1] , t );
         C_ADDTO( Fout[1] ,  t );

         t.r = Fout2[2].i;
         t.i = -Fout2[2].r;
         C_SUB( Fout2[2] ,  Fout[2] , t );
         C_ADDTO( Fout[2] ,  t );

         t.r = S_MUL(Fout2[3].i - Fout2[3].r, tw);
         t.i = S_MUL(-(Fout2[3].i + Fout2[3].r), tw);
         C_SUB( Fout2[3] ,  Fout[3] , t );
         C_ADDTO( Fout[3] ,  t );
         Fout += 8;
      }
   }
}

static void kf_bfly4(cfloat * Fout, const size_t fstride, const kiss_fft_state *st, int m, int N, int mm)
{
   int i;
   if (m==1){
      /* Degenerate case where all the twiddles are 1. */
      cfloat scratch0, scratch1;
      for (i=0;i<N;i++) {
         C_SUB( scratch0 , *Fout, Fout[2] );
         C_ADDTO(*Fout, Fout[2]);
         C_ADD( scratch1 , Fout[1] , Fout[3] );
         C_SUB( Fout[2], *Fout, scratch1 );
         C_ADDTO( *Fout , scratch1 );
         C_SUB( scratch1 , Fout[1] , Fout[3] );

         Fout[1].r = scratch0.r + scratch1.i;
         Fout[1].i = scratch0.i - scratch1.r;
         Fout[3].r = scratch0.r - scratch1.i;
         Fout[3].i = scratch0.i + scratch1.r;
         Fout+=4;
      }
   } else {
      int j;
      cfloat scratch[6];
      const cfloat *tw1,*tw2,*tw3;
      const int m2=2*m;
      const int m3=3*m;
      cfloat * Fout_beg = Fout;
      for (i=0;i<N;i++) {
         Fout = Fout_beg + i*mm;
         tw3 = tw2 = tw1 = st->twiddles;
         /* m is guaranteed to be a multiple of 4. */
         for (j=0;j<m;j++) {
            C_MUL(scratch[0],Fout[m] , *tw1 );
            C_MUL(scratch[1],Fout[m2] , *tw2 );
            C_MUL(scratch[2],Fout[m3] , *tw3 );

            C_SUB( scratch[5] , *Fout, scratch[1] );
            C_ADDTO(*Fout, scratch[1]);
            C_ADD( scratch[3] , scratch[0] , scratch[2] );
            C_SUB( scratch[4] , scratch[0] , scratch[2] );
            C_SUB( Fout[m2], *Fout, scratch[3] );
            tw1 += fstride;
            tw2 += fstride*2;
            tw3 += fstride*3;
            C_ADDTO( *Fout , scratch[3] );

            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
            ++Fout;
         }
      }
   }
}

#ifndef RADIX_TWO_ONLY
static void kf_bfly3(cfloat * Fout, const size_t fstride, const kiss_fft_state *st, int m, int N, int mm)
{
   int i;
   size_t k;
   const size_t m2 = 2*m;
   const cfloat *tw1,*tw2;
   cfloat scratch[5];
   cfloat epi3;

   cfloat * Fout_beg = Fout;
   epi3 = st->twiddles[fstride*m];
   for (i=0;i<N;i++) {
      Fout = Fout_beg + i*mm;
      tw1=tw2=st->twiddles;
      /* For non-custom modes, m is guaranteed to be a multiple of 4. */
      k=m;
      do {
         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;

         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);

         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;
         Fout[m].r = Fout[m].r - scratch[0].i;
         Fout[m].i = Fout[m].i + scratch[0].r;

         ++Fout;
      } while(--k);
   }
}

#ifndef OVERRIDE_kf_bfly5
static void kf_bfly5(cfloat * Fout, const size_t fstride, const kiss_fft_state *st, int m, int N, int mm)
{
   cfloat *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
   int i, u;
   cfloat scratch[13];
   const cfloat *tw;
   cfloat ya,yb;
   cfloat * Fout_beg = Fout;

   ya = st->twiddles[fstride*m];
   yb = st->twiddles[fstride*2*m];
   tw=st->twiddles;

   for (i=0;i<N;i++) {
      Fout = Fout_beg + i*mm;
      Fout0=Fout;
      Fout1=Fout0+m;
      Fout2=Fout0+2*m;
      Fout3=Fout0+3*m;
      Fout4=Fout0+4*m;

      /* For non-custom modes, m is guaranteed to be a multiple of 4. */
      for ( u=0; u<m; ++u ) {
         scratch[0] = *Fout0;

         C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
         C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
         C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
         C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

         C_ADD( scratch[7],scratch[1],scratch[4]);
         C_SUB( scratch[10],scratch[1],scratch[4]);
         C_ADD( scratch[8],scratch[2],scratch[3]);
         C_SUB( scratch[9],scratch[2],scratch[3]);

         Fout0->r = Fout0->r + (scratch[7].r + scratch[8].r);
         Fout0->i = Fout0->i + (scratch[7].i + scratch[8].i);

         scratch[5].r = scratch[0].r + (S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r));
         scratch[5].i = scratch[0].i + (S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r));

         scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
         scratch[6].i = -(S_MUL(scratch[10].r,ya.i) + S_MUL(scratch[9].r,yb.i));

         C_SUB(*Fout1,scratch[5],scratch[6]);
         C_ADD(*Fout4,scratch[5],scratch[6]);

         scratch[11].r = scratch[0].r + (S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r));
         scratch[11].i = scratch[0].i + (S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r));
         scratch[12].r = S_MUL(scratch[9].i,ya.i) - S_MUL(scratch[10].i,yb.i);
         scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

         C_ADD(*Fout2,scratch[11],scratch[12]);
         C_SUB(*Fout3,scratch[11],scratch[12]);

         ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
      }
   }
}
#endif /* OVERRIDE_kf_bfly5 */
#endif

#ifdef CUSTOM_MODES
static void compute_bitrev_table(int Fout, short *f, const size_t fstride, int in_stride, short * factors, const kiss_fft_state *st)
{
   const int p=*factors++; /* the radix  */
   const int m=*factors++; /* stage's fft length/p */
    /*printf ("fft %d %d %d %d %d %d\n", p*m, m, p, s2, fstride*in_stride, N);*/
   if (m==1) {
      int j;
      for (j=0;j<p;j++){
         *f = Fout+j;
         f += fstride*in_stride;
      }
   } else {
      int j;
      for (j=0;j<p;j++) {
         compute_bitrev_table( Fout , f, fstride*p, in_stride, factors,st);
         f += fstride*in_stride;
         Fout += m;
      }
   }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static int kf_factor(int n,short * facbuf)
{
    int p=4;
    int i;
    int stages=0;
    int nbak = n;

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p>32000 || (int)p*(int)p > n)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
#ifdef RADIX_TWO_ONLY
        if (p!=2 && p != 4)
#else
        if (p>5)
#endif
        {
           return 0;
        }
        facbuf[2*stages] = p;
        if (p==2 && stages > 1)
        {
           facbuf[2*stages] = 4;
           facbuf[2] = 2;
        }
        stages++;
    } while (n > 1);
    n = nbak;
    /* Reverse the order to get the radix 4 at the end, so we can use the
       fast degenerate case. It turns out that reversing the order also
       improves the noise behaviour. */
    for (i=0;i<stages/2;i++) {
       int tmp;
       tmp = facbuf[2*i];
       facbuf[2*i] = facbuf[2*(stages-i-1)];
       facbuf[2*(stages-i-1)] = tmp;
    }
    for (i=0;i<stages;i++) {
        n /= facbuf[2*i];
        facbuf[2*i+1] = n;
    }
    return 1;
}

static void compute_twiddles(cfloat *twiddles)
{
   int i;
   float phase;
   for (i=0;i<NFFT;++i) {
       phase = ( -2.0F*3.1415927F / NFFT) * i;
      kf_cexp(twiddles+i, phase );
   }

}

int opus_fft_alloc_arch_c(kiss_fft_state *st) {
   (void)st;
   return 0;
}

/*
 *
 * Allocates all necessary storage space for the fft and ifft.
 * The return value is a contiguous block of memory.  As such,
 * It can be freed with free().
 * */
kiss_fft_state KST;
cfloat TW[NFFT];
short BIT[NFFT];
kiss_fft_state *opus_fft_alloc_twiddles(void * mem,size_t * lenmem, const kiss_fft_state *base, int arch)
{
	kiss_fft_state *st= &KST;
	
    if ( lenmem!=NULL ) {
        if (mem != NULL && *lenmem >= sizeof(KST))
            st = (kiss_fft_state*)mem;
        *lenmem = sizeof(KST);
    }
    if (st) {
        short* bitrev = &BIT[0];
        cfloat *twiddles = &TW[0];

        st->nfft = NFFT;
        st->scale = 1.f/ NFFT;

        if (base != NULL)
        {
           st->twiddles = base->twiddles;
           st->shift = 0;
           while (st->shift < 32 && NFFT <<st->shift != base->nfft)
              st->shift++;
           if (st->shift>=32)
              return NULL;
        } else {
           st->twiddles = twiddles;
           compute_twiddles(twiddles);
           st->shift = -1;
        }
        if (!kf_factor(NFFT,st->factors))
           return NULL;
        /* bitrev */
        st->bitrev = bitrev;
        if (st->bitrev==NULL)
           return NULL;
        compute_bitrev_table(0, bitrev, 1,1, st->factors,st);

        /* Initialize architecture specific fft parameters */
        if (opus_fft_alloc_arch(st, arch))
           return NULL;
    }
    return st;
}

kiss_fft_state *opus_fft_alloc(void * mem,size_t * lenmem, int arch)
{
   return opus_fft_alloc_twiddles(mem, lenmem, NULL, arch);
}


#endif /* CUSTOM_MODES */

void opus_fft_impl(const kiss_fft_state *st,cfloat *fout)
{
    int m2, m;
    int p;
    int L;
    int fstride[MAXFACTORS];
    int i;
    int shift;

    /* st->shift can be -1 */
    shift = st->shift>0 ? st->shift : 0;

    fstride[0] = 1;
    L=0;
    do {
       p = st->factors[2*L];
       m = st->factors[2*L+1];
       fstride[L+1] = fstride[L]*p;
       L++;
    } while(m!=1);
    m = st->factors[2*L-1];
    for (i=L-1;i>=0;i--) {
       if (i!=0)
          m2 = st->factors[2*i-1];
       else
          m2 = 1;
       switch (st->factors[2*i]) {
       case 2:
          kf_bfly2(fout, m, fstride[i]);
          break;
       case 4:
          kf_bfly4(fout,fstride[i]<<shift,st,m, fstride[i], m2);
          break;
 #ifndef RADIX_TWO_ONLY
       case 3:
          kf_bfly3(fout,fstride[i]<<shift,st,m, fstride[i], m2);
          break;
       case 5:
          kf_bfly5(fout,fstride[i]<<shift,st,m, fstride[i], m2);
          break;
 #endif
       }
       m = m2;
    }
}

void opus_fft_c(const kiss_fft_state *st,const cfloat *fin,cfloat *fout)
{
   int i;
   float scale;

   scale = st->scale;

   /* Bit-reverse the input */
   for (i=0;i<st->nfft;i++) {
      cfloat x = fin[i];
      fout[st->bitrev[i]].r = scale * x.r;
      fout[st->bitrev[i]].i = scale * x.i;
   }
   opus_fft_impl(st, fout);
}


void opus_ifft_c(const kiss_fft_state *st,const cfloat *fin,cfloat *fout)
{
   int i;

   /* Bit-reverse the input */
   for (i=0;i<st->nfft;i++)
      fout[st->bitrev[i]] = fin[i];
   for (i=0;i<st->nfft;i++)
      fout[i].i = -fout[i].i;
   opus_fft_impl(st, fout);
   for (i=0;i<st->nfft;i++)
      fout[i].i = -fout[i].i;
}
