/*
    Copyright (C) 2008-2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"
#include "fft.h"
#include "fft_tuning.h"
#include "flint.h"

void _fmpz_poly_mullow_SS(fmpz * output, const fmpz * input1, slong len1, 
               const fmpz * input2, slong len2, slong trunc)
{
    slong len_out, loglen, loglen2, n;
    slong output_bits, limbs, size, i;
    mp_limb_t * ptr, ** t1, ** t2, ** tt, ** s1, ** ii, ** jj;
    slong bits1, bits2;
    ulong size1, size2;
    int sign = 0;
    int N;
    TMP_INIT;

    TMP_START;

    len1 = FLINT_MIN(len1, trunc);
    len2 = FLINT_MIN(len2, trunc);

    len_out = len1 + len2 - 1;
    loglen  = FLINT_CLOG2(len_out);
    loglen2 = FLINT_CLOG2(len2);
    n = (WORD(1) << (loglen - 2));

    size1 = _fmpz_vec_max_limbs(input1, len1); 
    size2 = _fmpz_vec_max_limbs(input2, len2);

    /* Start with an upper bound on the number of bits needed */
    output_bits = FLINT_BITS * (size1 + size2) + loglen2 + 1; 
    
    /* round up for sqrt2 trick */
    output_bits = (((output_bits - 1) >> (loglen - 2)) + 1) << (loglen - 2);

    limbs = (output_bits - 1) / FLINT_BITS + 1; /* initial size of FFT coeffs */
    if (limbs > FFT_MULMOD_2EXPP1_CUTOFF) /* can't be worse than next power of 2 limbs */
        limbs = (WORD(1) << FLINT_CLOG2(limbs));
    size = limbs + 1;

    /* allocate space for ffts */

    N = flint_get_num_threads();
    ii = flint_malloc((4*(n + n*size) + 5*size*N)*sizeof(mp_limb_t));
    for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
        ii[i] = ptr;
   t1 = TMP_ALLOC(N*sizeof(mp_limb_t *));
   t2 = TMP_ALLOC(N*sizeof(mp_limb_t *));
   s1 = TMP_ALLOC(N*sizeof(mp_limb_t *));
   tt = TMP_ALLOC(N*sizeof(mp_limb_t *));

   t1[0] = ptr;
   t2[0] = t1[0] + size*N;
   s1[0] = t2[0] + size*N;
   tt[0] = s1[0] + size*N;

   for (i = 1; i < N; i++)
   {
      t1[i] = t1[i - 1] + size;
      t2[i] = t2[i - 1] + size;
      s1[i] = s1[i - 1] + size;
      tt[i] = tt[i - 1] + 2*size;
   }

    if (input1 != input2)
    {
        jj = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
        for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
            jj[i] = ptr;
    } else jj = ii;

    /* put coefficients into FFT vecs */
    bits1 = _fmpz_vec_get_fft(ii, input1, limbs, len1);
    for (i = len1; i < 4*n; i++)
        flint_mpn_zero(ii[i], limbs + 1);

    if (input1 != input2) 
    {
        bits2 = _fmpz_vec_get_fft(jj, input2, limbs, len2);
        for (i = len2; i < 4*n; i++)
            flint_mpn_zero(jj[i], limbs + 1);
    }
    else bits2 = bits1;

    if (bits1 < WORD(0) || bits2 < WORD(0)) 
    {
        sign = 1;  
        bits1 = FLINT_ABS(bits1);
        bits2 = FLINT_ABS(bits2);
    }

    /* Recompute the number of bits/limbs now that we know how large everything is */
    output_bits = bits1 + bits2 + loglen2 + sign;

    /* round up output bits for sqrt2 */
    output_bits = (((output_bits - 1) >> (loglen - 2)) + 1) << (loglen - 2);

    limbs = (output_bits - 1) / FLINT_BITS + 1;
    limbs = fft_adjust_limbs(limbs); /* round up limbs for Nussbaumer */
    
    fft_convolution(ii, jj, loglen - 2, limbs, len_out, t1, t2, s1, tt); 

    _fmpz_vec_set_fft(output, trunc, ii, limbs, sign); /* write output */

    flint_free(ii); 
    if (input1 != input2) 
        flint_free(jj);

    TMP_END;
}

void
fmpz_poly_mullow_SS(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len1 <= 2 || len2 <= 2 || n <= 2)
    {
        fmpz_poly_mullow_classical(res, poly1, poly2, n);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);
    fmpz_poly_fit_length(res, n);

    if (len1 >= len2)
        _fmpz_poly_mullow_SS(res->coeffs, poly1->coeffs, len1,
                                          poly2->coeffs, len2, n);
    else
        _fmpz_poly_mullow_SS(res->coeffs, poly2->coeffs, len2,
                                          poly1->coeffs, len1, n);

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
