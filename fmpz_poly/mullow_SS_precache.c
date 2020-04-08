/*
    Copyright (C) 2008-2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"
#include "fft.h"
#include "fft_tuning.h"
#include "flint.h"

void fmpz_poly_mul_SS_precache_init(fmpz_poly_precache_t pre,
                 const fmpz_poly_t poly1, slong len2, slong bits2, slong trunc)
{
    slong i, len_out, loglen2;
    slong output_bits, size;
    ulong size1, size2;
    mp_limb_t * ptr;
    int N;

    pre->len1 = poly1->length;

    len_out = pre->len1 + len2 - 1;
    pre->loglen  = FLINT_CLOG2(len_out);
    loglen2 = FLINT_CLOG2(FLINT_MIN(pre->len1, len2));
    pre->n = (WORD(1) << (pre->loglen - 2));

    size1 = _fmpz_vec_max_limbs(poly1->coeffs, pre->len1);
    size2 = FLINT_ABS(bits2);
    size2 = (size2 + FLINT_BITS - 1)/FLINT_BITS;

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
    pre->ii = (mp_limb_t **)
        flint_malloc((4*(pre->n + pre->n*size) + 5*size*N)*sizeof(mp_limb_t));
    for (i = 0, ptr = (mp_limb_t *) pre->ii + 4*pre->n; i < 4*pre->n; i++, ptr += size) 
        pre->ii[i] = ptr;
    pre->t1 = (mp_limb_t **) flint_malloc(4*N*sizeof(mp_limb_t *));
    pre->t2 = (mp_limb_t **) pre->t1 + N;
    pre->s1 = (mp_limb_t **) pre->t2 + N;
    pre->tt = (mp_limb_t **) pre->s1 + N;

    pre->t1[0] = ptr;
    pre->t2[0] = pre->t1[0] + size*N;
    pre->s1[0] = pre->t2[0] + size*N;
    pre->tt[0] = pre->s1[0] + size*N;

    for (i = 1; i < N; i++)
    {
        pre->t1[i] = pre->t1[i - 1] + size;
        pre->t2[i] = pre->t2[i - 1] + size;
        pre->s1[i] = pre->s1[i - 1] + size;
        pre->tt[i] = pre->tt[i - 1] + 2*size;
    }

    /* put coefficients into FFT vecs */
    pre->bits1 = _fmpz_vec_get_fft(pre->ii, poly1->coeffs,
                                                        pre->limbs, pre->len1);
    for (i = pre->len1; i < 4*pre->n; i++)
        flint_mpn_zero(pre->ii[i], size);

    pre->bits1 = FLINT_ABS(pre->bits1);

    /* Recompute the number of bits/limbs now that we know how large everything is */
    output_bits = pre->bits1 + bits2 + loglen2 + 1;

    /* round up output bits for sqrt2 */
    output_bits = (((output_bits - 1) >> (pre->loglen - 2)) + 1) << (pre->loglen - 2);

    pre->limbs = (output_bits - 1) / FLINT_BITS + 1;
    pre->limbs = fft_adjust_limbs(pre->limbs); /* round up limbs for Nussbaumer */

    fft_precache(pre->ii, pre->loglen - 2, pre->limbs, trunc,
                                                    pre->t1, pre->t2, pre->s1);

    fmpz_poly_init(pre->poly1);
    fmpz_poly_set(pre->poly1, poly1);
}

void fmpz_poly_mul_precache_clear(fmpz_poly_precache_t pre)
{
    flint_free(pre->ii);
    flint_free(pre->t1);
    fmpz_poly_clear(pre->poly1);
}

void _fmpz_poly_mullow_SS_precache(fmpz * output, fmpz_poly_precache_t pre,
                                  const fmpz * input2, slong len2, slong trunc)
{
    slong len_out;
    slong size, i;
    mp_limb_t * ptr, ** jj;

    len_out = FLINT_MAX(pre->len1 + len2 - 1, 2*pre->n);

    size = pre->limbs + 1;

    /* allocate space for fft */

    jj = flint_malloc(4*(pre->n + pre->n*size)*sizeof(mp_limb_t));
    for (i = 0, ptr = (mp_limb_t *) jj + 4*pre->n; i < 4*pre->n; i++, ptr += size) 
        jj[i] = ptr;

    /* put coefficients into FFT vecs */
    _fmpz_vec_get_fft(jj, input2, pre->limbs, len2);
    for (i = len2; i < 4*pre->n; i++)
        flint_mpn_zero(jj[i], size);

    fft_convolution_precache(jj, pre->ii, pre->loglen - 2, pre->limbs,
                                  len_out, pre->t1, pre->t2, pre->s1, pre->tt);

    _fmpz_vec_set_fft(output, trunc, jj, pre->limbs, 1); /* write output */

    flint_free(jj);
}

void
fmpz_poly_mullow_SS_precache(fmpz_poly_t res,
                    fmpz_poly_precache_t pre, const fmpz_poly_t poly2, slong n)
{
    const slong len2 = poly2->length;

    if (pre->len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (pre->len1 <= 2 || len2 <= 2 || n <= 2)
    {
        fmpz_poly_mullow_classical(res, pre->poly1, poly2, n);
        return;
    }

    n = FLINT_MIN(n, pre->len1 + len2 - 1);
    fmpz_poly_fit_length(res, n);

    _fmpz_poly_mullow_SS_precache(res->coeffs, pre,
                                          poly2->coeffs, len2, n);

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
