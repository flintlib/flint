/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FFT_H
#define FFT_H

#ifdef FFT_INLINES_C
#define FFT_INLINE
#else
#define FFT_INLINE static inline
#endif

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

FFT_INLINE
mp_limb_t mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    mp_limb_t ret;
    mp_ptr t;

    if (n == 0)
        return 0;

    if ((s == x && d == y) || (s == y && d == x))
    {
        t = (mp_ptr) flint_malloc(n * sizeof(mp_limb_t));
        ret = mpn_sub_n(t, x, y, n);
        ret += 2 * mpn_add_n(s, x, y, n);
        flint_mpn_copyi(d, t, n);
        flint_free(t);
        return ret;
    }

    if (s == x || s == y)
    {
        ret = mpn_sub_n(d, x, y, n);
        ret += 2 * mpn_add_n(s, x, y, n);
        return ret;
    }

    ret = 2 * mpn_add_n(s, x, y, n);
    ret += mpn_sub_n(d, x, y, n);
    return ret;
}

#define fft_sumdiff(t, u, r, s, n) \
   (n == 0 ? 0 : mpn_sumdiff_n(t, u, r, s, n))


#define SWAP_PTRS(xx, yy) \
   do { \
      mp_limb_t * __ptr = xx; \
      xx = yy; \
      yy = __ptr; \
   } while (0)

/* used for generating random values mod p in test code */
/* NOTE: One needs to include mpn_extras.h for this macro. */
#define random_fermat(nn, state, limbs) \
   do { \
      if (n_randint(state, 10) == 0) { \
         flint_mpn_zero(nn, limbs); \
         nn[limbs] = 1; \
      } else { \
         if (n_randint(state, 2) == 0) \
            flint_mpn_rrandom(nn, state->gmp_state, limbs); \
         else \
            flint_mpn_urandomb(nn, state->gmp_state, limbs*FLINT_BITS); \
         nn[limbs] = n_randint(state, 1024); \
      } \
      if (n_randint(state, 2)) \
         nn[limbs] = -nn[limbs]; \
   } while (0)

FFT_INLINE
void mpn_addmod_2expp1_1(mp_limb_t * r, mp_size_t limbs, mp_limb_signed_t c)
{
   mp_limb_t sum = r[0] + c;

   /* check if adding c would cause a carry to propagate */
   if ((mp_limb_signed_t)(sum ^ r[0]) >= 0)
      r[0] = sum;
   else
   {
      if (c >= 0) mpn_add_1(r, r, limbs + 1, c);
      else mpn_sub_1(r, r, limbs + 1, -c);
   }
}

void fft_combine_limbs(mp_limb_t * res, mp_limb_t ** poly, slong length,
            mp_size_t coeff_limbs, mp_size_t output_limbs, mp_size_t total_limbs);

void fft_combine_bits(mp_limb_t * res, mp_limb_t ** poly, slong length,
                 flint_bitcnt_t bits, mp_size_t output_limbs, mp_size_t total_limbs);

mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_srcptr limbs,
            mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs);

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_srcptr limbs,
                 mp_size_t total_limbs, flint_bitcnt_t bits, mp_size_t output_limbs);

void fermat_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs);

void mpn_normmod_2expp1(mp_limb_t * t, mp_size_t limbs);

void butterfly_lshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1,
                       mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y);

void butterfly_rshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1,
                       mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y);

void mpn_negmod_2expp1(mp_limb_t* z, const mp_limb_t* a, mp_size_t limbs);

void mpn_mul_2expmod_2expp1(mp_limb_t * t,
                                  mp_limb_t * i1, mp_size_t limbs, flint_bitcnt_t d);

void mpn_div_2expmod_2expp1(mp_limb_t * t,
                                  mp_limb_t * i1, mp_size_t limbs, flint_bitcnt_t d);

void fft_adjust(mp_limb_t * r, mp_limb_t * i1,
                                     mp_size_t i, mp_size_t limbs, flint_bitcnt_t w);

void fft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
                     mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w);

void ifft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
                     mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w);

void fft_radix2(mp_limb_t ** ii,
                    mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2);

void fft_truncate1(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
                               mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void fft_truncate(mp_limb_t ** ii,  mp_size_t n, flint_bitcnt_t w,
                               mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void ifft_radix2(mp_limb_t ** ii, mp_size_t n,
                                 flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2);

void ifft_truncate1(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
                               mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void ifft_truncate(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
                               mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void fft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t,
                         mp_limb_t * i1, mp_limb_t * i2, mp_size_t i,
                                mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp);

void ifft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
   mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp);

void fft_adjust_sqrt2(mp_limb_t * r, mp_limb_t * i1,
                   mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp);

void fft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
            mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc);

void ifft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
            mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc);

void mul_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
                        mp_srcptr i2, mp_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w);

void fft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
   mp_limb_t * s, mp_limb_t * t, mp_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2);

void ifft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
   mp_limb_t * s, mp_limb_t * t, mp_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2);

void fft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
      mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                            mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs);

void ifft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
        mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                            mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs);

void fft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
        mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
           mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc);

void ifft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
        mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
           mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc);

void fft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n,
                       flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                                mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc);

void ifft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n,
                      flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                                mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc);

void mul_mfa_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
                        mp_srcptr i2, mp_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w);

void fft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n,
                      flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                                mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc);

void fft_mfa_truncate_sqrt2_inner(mp_limb_t ** ii, mp_limb_t ** jj,
            mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc, mp_limb_t ** tt);

void ifft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n,
                        flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
                                mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc);

void fft_negacyclic(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
                             mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp);

void ifft_negacyclic(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
                             mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp);

void fft_naive_convolution_1(mp_limb_t * r, mp_limb_t * ii,
                                                     mp_limb_t * jj, mp_size_t m);

void _fft_mulmod_2expp1(mp_limb_t * r1, mp_limb_t * i1, mp_limb_t * i2,
                             mp_size_t r_limbs, flint_bitcnt_t depth, flint_bitcnt_t w);

slong fft_adjust_limbs(mp_size_t limbs);

void fft_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2,
                                        mp_size_t n, mp_size_t w, mp_limb_t * tt);

void flint_mpn_mul_fft_main(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
                        mp_srcptr i2, mp_size_t n2);

void fft_convolution_basic(mp_limb_t ** ii, mp_limb_t ** jj,
		     slong depth, slong limbs, slong trunc, mp_limb_t ** t1,
                            mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt);

void fft_convolution(mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
                                    slong limbs, slong trunc, mp_limb_t ** t1,
                            mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt);

/***** FFT Precaching *****/

void fft_precache(mp_limb_t ** jj, slong depth, slong limbs,
               slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1);

void fft_convolution_precache(mp_limb_t ** ii, mp_limb_t ** jj,
               slong depth, slong limbs, slong trunc, mp_limb_t ** t1,
	                    mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt);

#ifdef __cplusplus
}
#endif

#endif

