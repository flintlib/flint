/*
    Copyright (C) 2008-2011, 2020 William Hart

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

void fmpz_poly_mul_SS_precache_init(fmpz_poly_mul_precache_t pre,
                              slong len1, slong bits1, const fmpz_poly_t poly2)
{
    slong i, len_out, loglen2;
    slong output_bits, size;
    ulong size1, size2;
    mp_limb_t * ptr;
    mp_limb_t ** t1, ** t2, ** s1;
    int N;

    pre->len2 = poly2->length;

    len_out = len1 + pre->len2 - 1;
    pre->loglen  = FLINT_CLOG2(len_out);
    loglen2 = FLINT_CLOG2(FLINT_MIN(len1, pre->len2));
    pre->n = (WORD(1) << (pre->loglen - 2));
    size1 = FLINT_ABS(bits1);
    size1 = (size1 + FLINT_BITS - 1)/FLINT_BITS;
    size2 = _fmpz_vec_max_limbs(poly2->coeffs, pre->len2);

    /* Start with an upper bound on the number of bits needed */
    output_bits = FLINT_BITS*(size1 + size2) + loglen2 + 1; 

    /* round up for sqrt2 trick */
    output_bits = (((output_bits - 1) >> (pre->loglen - 2)) + 1) << (pre->loglen - 2);

    pre->limbs = (output_bits - 1) / FLINT_BITS + 1; /* initial size of FFT coeffs */

    if (pre->limbs > FFT_MULMOD_2EXPP1_CUTOFF) /* can't be worse than next power of 2 limbs */
        pre->limbs = (WORD(1) << FLINT_CLOG2(pre->limbs));
    size = pre->limbs + 1;

    /* allocate space for ffts */
    N = flint_get_num_threads();
    pre->jj = (mp_limb_t **)
        flint_malloc((4*(pre->n + pre->n*size) + 3*size*N + 3*N)*sizeof(mp_limb_t));
    for (i = 0, ptr = (mp_limb_t *) pre->jj + 4*pre->n; i < 4*pre->n; i++, ptr += size) 
        pre->jj[i] = ptr;
    t1 = (mp_limb_t **) ptr;
    t2 = (mp_limb_t **) t1 + N;
    s1 = (mp_limb_t **) t2 + N;
    ptr += 3*N;

    t1[0] = ptr;
    t2[0] = t1[0] + size*N;
    s1[0] = t2[0] + size*N;

    for (i = 1; i < N; i++)
    {
        t1[i] = t1[i - 1] + size;
        t2[i] = t2[i - 1] + size;
        s1[i] = s1[i - 1] + size;
    }

    /* put coefficients into FFT vecs */
    pre->bits2 = _fmpz_vec_get_fft(pre->jj, poly2->coeffs,
                                                        pre->limbs, pre->len2);
    for (i = pre->len2; i < 4*pre->n; i++)
        flint_mpn_zero(pre->jj[i], size);

    pre->bits2 = FLINT_ABS(pre->bits2);

    /* Recompute the number of bits/limbs now that we know how large everything is */
    output_bits = bits1 + pre->bits2 + loglen2 + 1;

    /* round up output bits for sqrt2 */
    output_bits = (((output_bits - 1) >> (pre->loglen - 2)) + 1) << (pre->loglen - 2);

    pre->limbs = (output_bits - 1) / FLINT_BITS + 1;
    pre->limbs = fft_adjust_limbs(pre->limbs); /* round up limbs for Nussbaumer */

    fft_precache(pre->jj, pre->loglen - 2, pre->limbs, len_out, t1, t2, s1);

    fmpz_poly_init(pre->poly2);
    fmpz_poly_set(pre->poly2, poly2);
}

void fmpz_poly_mul_precache_clear(fmpz_poly_mul_precache_t pre)
{
    flint_free(pre->jj);
    fmpz_poly_clear(pre->poly2);
}

void _fmpz_poly_mullow_SS_precache(fmpz * output, const fmpz * input1,
                         slong len1, fmpz_poly_mul_precache_t pre, slong trunc)
{
    slong len_out;
    slong size, i;
    mp_limb_t ** ii, ** t1, ** t2, ** s1, ** tt;
    mp_limb_t * ptr;
    int N;

    len_out = FLINT_MAX(len1 + pre->len2 - 1, 2*pre->n + 1);

    size = pre->limbs + 1;

    /* allocate space for ffts */
    N = flint_get_num_threads();
    ii = (mp_limb_t **)
        flint_malloc((4*(pre->n + pre->n*size) + 5*size*N + 4*N)*sizeof(mp_limb_t));
    for (i = 0, ptr = (mp_limb_t *) ii + 4*pre->n; i < 4*pre->n; i++, ptr += size) 
        ii[i] = ptr;
    t1 = (mp_limb_t **) ptr;
    t2 = (mp_limb_t **) t1 + N;
    s1 = (mp_limb_t **) t2 + N;
    tt = (mp_limb_t **) s1 + N;
    ptr += 4*N;

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

    /* put coefficients into FFT vecs */
    _fmpz_vec_get_fft(ii, input1, pre->limbs, len1);
    for (i = len1; i < 4*pre->n; i++)
        flint_mpn_zero(ii[i], size);

    fft_convolution_precache(ii, pre->jj, pre->loglen - 2, pre->limbs,
                                                      len_out, t1, t2, s1, tt);

    _fmpz_vec_set_fft(output, trunc, ii, pre->limbs, 1); /* write output */
    flint_free(ii);
}

void
fmpz_poly_mullow_SS_precache(fmpz_poly_t res,
                    const fmpz_poly_t poly1, fmpz_poly_mul_precache_t pre, slong n)
{
    const slong len1 = poly1->length;

    if (pre->len2 == 0 || len1 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (pre->len2 <= 2 || len1 <= 2 || n <= 2)
    {
        fmpz_poly_mullow_classical(res, poly1, pre->poly2, n);
        return;
    }

    n = FLINT_MIN(n, len1 + pre->len2 - 1);
    fmpz_poly_fit_length(res, n);

    _fmpz_poly_mullow_SS_precache(res->coeffs, poly1->coeffs, len1, pre, n);

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
