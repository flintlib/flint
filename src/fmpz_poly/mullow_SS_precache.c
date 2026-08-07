/*
    Copyright (C) 2008-2011, 2020 William Hart
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "mpn_extras.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft.h"
#include "fft_small.h"

void fmpz_poly_mul_SS_precache_init(fmpz_poly_mul_precache_t pre,
                              slong len1, slong bits1, const fmpz_poly_t poly2)
{
    slong negm;
    slong i, len_out, loglen2;
    slong output_bits, size;
    ulong * ptr;
    ulong ** t1, ** t2, ** s1;
    int N;

    pre->len2 = poly2->length;
    pre->bits2 = _fmpz_vec_max_bits(poly2->coeffs, pre->len2);
    pre->bits2 = FLINT_ABS(pre->bits2);

    len_out = len1 + pre->len2 - 1;
    pre->loglen  = FLINT_CLOG2(len_out);
    loglen2 = FLINT_CLOG2(FLINT_MIN(len1, pre->len2));
    pre->n = (WORD(1) << (pre->loglen - 2));

    /* The transform ring is fixed once, up front: the output bit
       bound uses the actual bit counts (both available here), and
       ring selection happens before any allocation so that the
       buffers, the packed operand, the transforms and the optional
       sd_fft pointwise context all live in the same ring. */
    output_bits = bits1 + pre->bits2 + loglen2 + 1;

    /* round up for sqrt2 trick */
    output_bits = (((output_bits - 1) >> (pre->loglen - 2)) + 1) << (pre->loglen - 2);

    pre->limbs = (output_bits - 1) / FLINT_BITS + 1;
    pre->limbs = fft_adjust_limbs_sd_fft(pre->limbs, pre->n, &negm);
    pre->sdctx = NULL;
#if FLINT_HAVE_FFT_SMALL
    if (negm != 0)
    {
        pre->sdctx = flint_malloc(sizeof(sd_fft_mpn_mulmod_2expp1_ctx_struct));
        sd_fft_mpn_mulmod_2expp1_ctx_init(pre->sdctx, get_default_mpn_ctx(),
                            FLINT_BITS*pre->limbs, negm);
    }
#else
    (void) negm;
#endif
    size = pre->limbs + 1;

    /* allocate space for ffts */
    N = flint_get_num_threads();
    pre->jj = (ulong **)
        flint_malloc((4*(pre->n + pre->n*size) + 3*size*N + 3*N)*sizeof(ulong));
    for (i = 0, ptr = (ulong *) pre->jj + 4*pre->n; i < 4*pre->n; i++, ptr += size)
        pre->jj[i] = ptr;
    t1 = (ulong **) ptr;
    t2 = (ulong **) t1 + N;
    s1 = (ulong **) t2 + N;
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
    _fmpz_vec_get_fft(pre->jj, poly2->coeffs, pre->limbs, pre->len2);
    for (i = pre->len2; i < 4*pre->n; i++)
        flint_mpn_zero(pre->jj[i], size);

    fft_precache(pre->jj, pre->loglen - 2, pre->limbs, len_out, t1, t2, s1);

    fmpz_poly_init(pre->poly2);
    fmpz_poly_set(pre->poly2, poly2);
}

void fmpz_poly_mul_precache_clear(fmpz_poly_mul_precache_t pre)
{
#if FLINT_HAVE_FFT_SMALL
    if (pre->sdctx != NULL)
    {
        sd_fft_mpn_mulmod_2expp1_ctx_clear(pre->sdctx);
        flint_free(pre->sdctx);
        pre->sdctx = NULL;
    }
#endif
    flint_free(pre->jj);
    fmpz_poly_clear(pre->poly2);
}

void _fmpz_poly_mullow_SS_precache(fmpz * output, const fmpz * input1,
                         slong len1, fmpz_poly_mul_precache_t pre, slong trunc)
{
    slong len_out;
    slong size, i;
    ulong ** ii, ** t1, ** t2, ** s1, ** tt;
    ulong * ptr;
    int N;

    len_out = FLINT_MAX(len1 + pre->len2 - 1, 2*pre->n + 1);

    size = pre->limbs + 1;

    /* allocate space for ffts */
    N = flint_get_num_threads();
    ii = (ulong **)
        flint_malloc((4*(pre->n + pre->n*size) + 5*size*N + 4*N)*sizeof(ulong));
    for (i = 0, ptr = (ulong *) ii + 4*pre->n; i < 4*pre->n; i++, ptr += size)
        ii[i] = ptr;
    t1 = (ulong **) ptr;
    t2 = (ulong **) t1 + N;
    s1 = (ulong **) t2 + N;
    tt = (ulong **) s1 + N;
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

    fft_convolution_precache_sd_fft(pre->sdctx, ii, pre->jj, pre->loglen - 2, pre->limbs,
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
