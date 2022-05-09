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
#define FFT_INLINE FLINT_DLL
#else
#define FFT_INLINE static __inline__
#endif

#include "flint.h"

#if HAVE_OPENMP
#include <omp.h> /* must come after flint.h */
#endif

#ifdef __cplusplus
 extern "C" {
#endif

#define fft_sumdiff(t, u, r, s, n) \
   (n == 0 ? 0 : mpn_sumdiff_n(t, u, r, s, n))

#define SWAP_PTRS(xx, yy) \
   do { \
      ulong * __ptr = xx; \
      xx = yy; \
      yy = __ptr; \
   } while (0)

/* used for generating random values mod p in test code */
#define random_fermat(nn, state, limbs) \
   do { \
      if (n_randint(state, 10) == 0) { \
         FLINT_MPN_ZERO(nn, limbs); \
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

FLINT_DLL void fft_combine_limbs(ulong * res, ulong ** poly, slong length, 
            mp_mock_size_t coeff_limbs, mp_mock_size_t output_limbs, mp_mock_size_t total_limbs);

FLINT_DLL void fft_combine_bits(ulong * res, ulong ** poly, slong length, 
                 flint_bitcnt_t bits, mp_mock_size_t output_limbs, mp_mock_size_t total_limbs);

FLINT_DLL mp_mock_size_t fft_split_limbs(ulong ** poly, ulong_srcptr limbs, 
            mp_mock_size_t total_limbs, mp_mock_size_t coeff_limbs, mp_mock_size_t output_limbs);

FLINT_DLL mp_mock_size_t fft_split_bits(ulong ** poly, ulong_srcptr limbs, 
                 mp_mock_size_t total_limbs, flint_bitcnt_t bits, mp_mock_size_t output_limbs);

#ifdef __GMP_H__
FLINT_DLL void fermat_to_mpz(mpz_t m, ulong * i, mp_mock_size_t limbs);
#endif

FLINT_DLL void mpn_normmod_2expp1(ulong * t, mp_mock_size_t limbs);

FLINT_DLL void butterfly_lshB(ulong * t, ulong * u, ulong * i1, 
                       ulong * i2, mp_mock_size_t limbs, mp_mock_size_t x, mp_mock_size_t y);

FLINT_DLL void butterfly_rshB(ulong * t, ulong * u, ulong * i1, 
                       ulong * i2, mp_mock_size_t limbs, mp_mock_size_t x, mp_mock_size_t y);

FLINT_DLL void mpn_mul_2expmod_2expp1(ulong * t, 
                                  ulong * i1, mp_mock_size_t limbs, flint_bitcnt_t d);

FLINT_DLL void mpn_div_2expmod_2expp1(ulong * t, 
                                  ulong * i1, mp_mock_size_t limbs, flint_bitcnt_t d);

FLINT_DLL void fft_adjust(ulong * r, ulong * i1, 
                                     mp_mock_size_t i, mp_mock_size_t limbs, flint_bitcnt_t w);

FLINT_DLL void fft_butterfly(ulong * s, ulong * t, ulong * i1, 
                     ulong * i2, mp_mock_size_t i, mp_mock_size_t limbs, flint_bitcnt_t w);

FLINT_DLL void ifft_butterfly(ulong * s, ulong * t, ulong * i1, 
                     ulong * i2, mp_mock_size_t i, mp_mock_size_t limbs, flint_bitcnt_t w);

FLINT_DLL void fft_radix2(ulong ** ii, 
                    mp_mock_size_t n, flint_bitcnt_t w, ulong ** t1, ulong ** t2);

FLINT_DLL void fft_truncate1(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
                               ulong ** t1, ulong ** t2, mp_mock_size_t trunc);

FLINT_DLL void fft_truncate(ulong ** ii,  mp_mock_size_t n, flint_bitcnt_t w, 
                               ulong ** t1, ulong ** t2, mp_mock_size_t trunc);

FLINT_DLL void ifft_radix2(ulong ** ii, mp_mock_size_t n, 
                                 flint_bitcnt_t w, ulong ** t1, ulong ** t2);

FLINT_DLL void ifft_truncate1(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
                               ulong ** t1, ulong ** t2, mp_mock_size_t trunc);

FLINT_DLL void ifft_truncate(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
                               ulong ** t1, ulong ** t2, mp_mock_size_t trunc);

FLINT_DLL void fft_butterfly_sqrt2(ulong * s, ulong * t, 
                         ulong * i1, ulong * i2, mp_mock_size_t i, 
                                mp_mock_size_t limbs, flint_bitcnt_t w, ulong * temp);

FLINT_DLL void ifft_butterfly_sqrt2(ulong * s, ulong * t, ulong * i1, 
   ulong * i2, mp_mock_size_t i, mp_mock_size_t limbs, flint_bitcnt_t w, ulong * temp);

FLINT_DLL void fft_adjust_sqrt2(ulong * r, ulong * i1, 
                   mp_mock_size_t i, mp_mock_size_t limbs, flint_bitcnt_t w, ulong * temp);

FLINT_DLL void fft_truncate_sqrt2(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
            ulong ** t1, ulong ** t2, ulong ** temp, mp_mock_size_t trunc);

FLINT_DLL void ifft_truncate_sqrt2(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
            ulong ** t1, ulong ** t2, ulong ** temp, mp_mock_size_t trunc);

FLINT_DLL void mul_truncate_sqrt2(ulong_ptr r1, ulong_srcptr i1, mp_mock_size_t n1,
                        ulong_srcptr i2, mp_mock_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w);

FLINT_DLL void fft_butterfly_twiddle(ulong * u, ulong * v, 
   ulong * s, ulong * t, mp_mock_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2);

FLINT_DLL void ifft_butterfly_twiddle(ulong * u, ulong * v, 
   ulong * s, ulong * t, mp_mock_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2);

FLINT_DLL void fft_radix2_twiddle(ulong ** ii, mp_mock_size_t is,
      mp_mock_size_t n, flint_bitcnt_t w, ulong ** t1, ulong ** t2,
                            mp_mock_size_t ws, mp_mock_size_t r, mp_mock_size_t c, mp_mock_size_t rs);

FLINT_DLL void ifft_radix2_twiddle(ulong ** ii, mp_mock_size_t is,
        mp_mock_size_t n, flint_bitcnt_t w, ulong ** t1, ulong ** t2,
                            mp_mock_size_t ws, mp_mock_size_t r, mp_mock_size_t c, mp_mock_size_t rs);

FLINT_DLL void fft_truncate1_twiddle(ulong ** ii, mp_mock_size_t is,
        mp_mock_size_t n, flint_bitcnt_t w, ulong ** t1, ulong ** t2,
           mp_mock_size_t ws, mp_mock_size_t r, mp_mock_size_t c, mp_mock_size_t rs, mp_mock_size_t trunc);

FLINT_DLL void ifft_truncate1_twiddle(ulong ** ii, mp_mock_size_t is,
        mp_mock_size_t n, flint_bitcnt_t w, ulong ** t1, ulong ** t2,
           mp_mock_size_t ws, mp_mock_size_t r, mp_mock_size_t c, mp_mock_size_t rs, mp_mock_size_t trunc);

FLINT_DLL void fft_mfa_truncate_sqrt2(ulong ** ii, mp_mock_size_t n, 
                       flint_bitcnt_t w, ulong ** t1, ulong ** t2, 
                                ulong ** temp, mp_mock_size_t n1, mp_mock_size_t trunc);

FLINT_DLL void ifft_mfa_truncate_sqrt2(ulong ** ii, mp_mock_size_t n, 
                      flint_bitcnt_t w, ulong ** t1, ulong ** t2,
                                ulong ** temp, mp_mock_size_t n1, mp_mock_size_t trunc);

FLINT_DLL void mul_mfa_truncate_sqrt2(ulong_ptr r1, ulong_srcptr i1, mp_mock_size_t n1,
                        ulong_srcptr i2, mp_mock_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w);

FLINT_DLL void fft_mfa_truncate_sqrt2_outer(ulong ** ii, mp_mock_size_t n, 
                      flint_bitcnt_t w, ulong ** t1, ulong ** t2, 
                                ulong ** temp, mp_mock_size_t n1, mp_mock_size_t trunc);

FLINT_DLL void fft_mfa_truncate_sqrt2_inner(ulong ** ii, ulong ** jj, 
            mp_mock_size_t n, flint_bitcnt_t w, ulong ** t1, ulong ** t2, 
                ulong ** temp, mp_mock_size_t n1, mp_mock_size_t trunc, ulong ** tt);

FLINT_DLL void ifft_mfa_truncate_sqrt2_outer(ulong ** ii, mp_mock_size_t n, 
                        flint_bitcnt_t w, ulong ** t1, ulong ** t2, 
                                ulong ** temp, mp_mock_size_t n1, mp_mock_size_t trunc);

FLINT_DLL void fft_negacyclic(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
                             ulong ** t1, ulong ** t2, ulong ** temp);

FLINT_DLL void ifft_negacyclic(ulong ** ii, mp_mock_size_t n, flint_bitcnt_t w, 
                             ulong ** t1, ulong ** t2, ulong ** temp);

FLINT_DLL void fft_naive_convolution_1(ulong * r, ulong * ii, 
                                                     ulong * jj, mp_mock_size_t m);

FLINT_DLL void _fft_mulmod_2expp1(ulong * r1, ulong * i1, ulong * i2, 
                             mp_mock_size_t r_limbs, flint_bitcnt_t depth, flint_bitcnt_t w);

FLINT_DLL slong fft_adjust_limbs(mp_mock_size_t limbs);

FLINT_DLL void fft_mulmod_2expp1(ulong * r, ulong * i1, ulong * i2, 
                                        mp_mock_size_t n, mp_mock_size_t w, ulong * tt);

FLINT_DLL void flint_mpn_mul_fft_main(ulong_ptr r1, ulong_srcptr i1, mp_mock_size_t n1,
                        ulong_srcptr i2, mp_mock_size_t n2);

FLINT_DLL void fft_convolution_basic(ulong ** ii, ulong ** jj,
		     slong depth, slong limbs, slong trunc, ulong ** t1, 
                            ulong ** t2, ulong ** s1, ulong ** tt);
	
FLINT_DLL void fft_convolution(ulong ** ii, ulong ** jj, slong depth, 
                                    slong limbs, slong trunc, ulong ** t1, 
                            ulong ** t2, ulong ** s1, ulong ** tt);

/***** FFT Precaching *****/

FLINT_DLL void fft_precache(ulong ** jj, slong depth, slong limbs,
               slong trunc, ulong ** t1, ulong ** t2, ulong ** s1);

FLINT_DLL void fft_convolution_precache(ulong ** ii, ulong ** jj,
               slong depth, slong limbs, slong trunc, ulong ** t1,
	                    ulong ** t2, ulong ** s1, ulong ** tt);

#ifdef __cplusplus
}
#endif

#endif
