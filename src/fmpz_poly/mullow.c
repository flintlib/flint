/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

static void
_fmpz_poly_mullow_tiny1(fmpz * res, const fmpz * poly1,
                         slong len1, const fmpz * poly2, slong len2, slong n)
{
    slong i, j, c;

    _fmpz_vec_zero(res, n);

    for (i = 0; i < len1; i++)
    {
        c = poly1[i];

        if (c != 0)
        {
            for (j = 0; j < FLINT_MIN(len2, n - i); j++)
                res[i + j] += c * poly2[j];
        }
    }
}

static void
_fmpz_poly_mullow_tiny2(fmpz * res, const fmpz * poly1,
                         slong len1, const fmpz * poly2, slong len2, slong n)
{
    slong i, j, k, c, d;
    ulong hi, lo;
    nn_ptr tmp;
    TMP_INIT;

    TMP_START;

    tmp = TMP_ALLOC(2 * n * sizeof(ulong));

    flint_mpn_zero(tmp, 2 * n);

    for (i = 0; i < len1; i++)
    {
        c = poly1[i];

        if (c != 0)
        {
            for (j = 0; j < FLINT_MIN(len2, n - i); j++)
            {
                k = i + j;

                d = poly2[j];

                if (d != 0)
                {
                    smul_ppmm(hi, lo, c, d);
                    add_ssaaaa(tmp[2 * k + 1], tmp[2 * k],
                               tmp[2 * k + 1], tmp[2 * k], hi, lo);
                }
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        lo = tmp[2 * i];
        hi = tmp[2 * i + 1];

        if (((slong) hi) >= 0)
        {
            fmpz_set_uiui(res + i, hi, lo);
        }
        else
        {
            sub_ddmmss(hi, lo, 0, 0, hi, lo);
            fmpz_neg_uiui(res + i, hi, lo);
        }
    }

    TMP_END;
    return;
}

void
_fmpz_poly_mullow(fmpz * res, const fmpz * poly1, slong len1,
                                const fmpz * poly2, slong len2, slong n)
{
    slong bits1, bits2, rbits;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(len1 > 0);
    FLINT_ASSERT(len2 > 0);
    FLINT_ASSERT(n <= len1 + len2 - 1);

    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1, len1, poly2);
        return;
    }

    if (len1 < len2)
    {
        _fmpz_poly_mullow(res, poly2, len2, poly1, len1, n);
        return;
    }

    if (poly1 == poly2 && len1 == len2)
    {
        _fmpz_poly_sqrlow(res, poly1, len1, n);
        return;
    }

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits1 = FLINT_ABS(bits1);
    bits2 = FLINT_ABS(bits2);

#if FLINT_HAVE_FFT_SMALL
    if (len2 >= 100 && (bits1 + bits2 <= 40 || bits1 + bits2 >= 128 || len2 >= 200))
        if (_fmpz_poly_mul_mid_default_mpn_ctx(res, 0, n, poly1, len1, poly2, len2))
            return;
#endif

    if (bits1 <= SMALL_FMPZ_BITCOUNT_MAX && bits2 <= SMALL_FMPZ_BITCOUNT_MAX &&
        (len2 < 50 || (4 * len2 >= 3 * n && n < 150 + bits1 + bits2)))
    {
        rbits = bits1 + bits2 + FLINT_BIT_COUNT(len2);

        if (rbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_poly_mullow_tiny1(res, poly1, len1, poly2, len2, n);
            return;
        }
        else if (rbits <= 2 * FLINT_BITS - 1)
        {
            _fmpz_poly_mullow_tiny2(res, poly1, len1, poly2, len2, n);
            return;
        }
    }

#if FLINT_HAVE_FFT_SMALL

    /* same as in mul.c */
    if (len2 <= 6 && FLINT_MIN(bits1, bits2) <= 5000)
        _fmpz_poly_mullow_classical(res, poly1, len1, poly2, len2, n);
    else if (len2 <= 4 || (len2 <= 8 && bits1 + bits2 >= 1500 && bits1 + bits2 <= 10000))
        _fmpz_poly_mullow_karatsuba(res, poly1, len1, poly2, len2, n);
    else if
        ((len2 >= 8 && len2 <= 75 && bits1 + bits2 >= 800 && bits1 + bits2 <= 4000) ||
            (len1 + len2 >= 5000 && bits1 + bits2 >= 5000 + (len1 + len2) / 10 && flint_get_num_threads() >= 4))
        _fmpz_poly_mullow_SS(res, poly1, len1, poly2, len2, n);
    else
        _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);

#else

    if (len2 < 7)
    {
        _fmpz_poly_mullow_classical(res, poly1, len1, poly2, len2, n);
    }
    else
    {
        slong limbs1, limbs2;

        limbs1 = (bits1 + FLINT_BITS - 1) / FLINT_BITS;
        limbs2 = (bits2 + FLINT_BITS - 1) / FLINT_BITS;

        if (n < 16 && (limbs1 > 12 || limbs2 > 12))
            _fmpz_poly_mullow_karatsuba(res, poly1, len1, poly2, len2, n);
        else if (limbs1 + limbs2 <= 8)
            _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);
        else if ((limbs1+limbs2)/2048 > len1 + len2)
            _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);
        else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
            _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);
        else
            _fmpz_poly_mullow_SS(res, poly1, len1, poly2, len2, n);
    }

#endif
}

void
fmpz_poly_mullow(fmpz_poly_t res,
                   const fmpz_poly_t poly1, const fmpz_poly_t poly2,
                   slong n)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_mullow(t, poly1, poly2, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);

    fmpz_poly_fit_length(res, n);
    if (len1 >= len2)
        _fmpz_poly_mullow(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, n);
    else
        _fmpz_poly_mullow(res->coeffs, poly2->coeffs, len2, poly1->coeffs, len1, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
