/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2014 Fredrik Johansson

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
# include "fft_small.h"
#endif

static void
_fmpz_poly_mul_tiny1(fmpz * res, const fmpz * poly1,
                         slong len1, const fmpz * poly2, slong len2)
{
    slong i, j, c;

    _fmpz_vec_zero(res, len1 + len2 - 1);

    for (i = 0; i < len1; i++)
    {
        c = poly1[i];

        if (c != 0)
        {
            for (j = 0; j < len2; j++)
                res[i + j] += c * poly2[j];
        }
    }
}

static void
_fmpz_poly_mul_tiny2(fmpz * res, const fmpz * poly1,
                         slong len1, const fmpz * poly2, slong len2)
{
    slong i, j, k, c, d;
    mp_limb_t hi, lo;
    mp_ptr tmp;
    TMP_INIT;

    TMP_START;

    tmp = TMP_ALLOC(2 * (len1 + len2 - 1) * sizeof(mp_limb_t));

    flint_mpn_zero(tmp, 2 * (len1 + len2 - 1));

    for (i = 0; i < len1; i++)
    {
        c = poly1[i];

        if (c != 0)
        {
            for (j = 0; j < len2; j++)
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

    for (i = 0; i < len1 + len2 - 1; i++)
    {
        lo = tmp[2 * i];
        hi = tmp[2 * i + 1];

        if (((mp_limb_signed_t) hi) >= 0)
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
_fmpz_poly_mul(fmpz * res, const fmpz * poly1,
                         slong len1, const fmpz * poly2, slong len2)
{
    slong bits1, bits2, rbits;

    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1, len1, poly2);
        return;
    }

    if (poly1 == poly2 && len1 == len2)
    {
        _fmpz_poly_sqr(res, poly1, len1);
        return;
    }

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits1 = FLINT_ABS(bits1);
    bits2 = FLINT_ABS(bits2);

#if FLINT_HAVE_FFT_SMALL
    if (len2 >= 80 && (bits1 + bits2 <= 40 || bits1 + bits2 >= 128 || len2 >= 100))
        if (_fmpz_poly_mul_mid_default_mpn_ctx(res, 0, len1 + len2 - 1, poly1, len1, poly2, len2))
            return;
#endif

    if (bits1 <= SMALL_FMPZ_BITCOUNT_MAX && bits2 <= SMALL_FMPZ_BITCOUNT_MAX &&
        (len2 < 40 + (bits1 + bits2) / 2 || len1 < 70 + (bits1 + bits2) / 2))
    {
        rbits = bits1 + bits2 + FLINT_BIT_COUNT(len2);

        if (rbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_poly_mul_tiny1(res, poly1, len1, poly2, len2);
            return;
        }
        else if (rbits <= 2 * FLINT_BITS - 1)
        {
            _fmpz_poly_mul_tiny2(res, poly1, len1, poly2, len2);
            return;
        }
    }

#if FLINT_HAVE_FFT_SMALL

    if (len2 <= 6 && FLINT_MIN(bits1, bits2) <= 5000)
        _fmpz_poly_mul_classical(res, poly1, len1, poly2, len2);
    else if (len2 <= 4 || (len2 <= 8 && bits1 + bits2 >= 1500 && bits1 + bits2 <= 10000))
        _fmpz_poly_mul_karatsuba(res, poly1, len1, poly2, len2);
    else if
        /* The fft_small-based KS is so efficient that SS currently
           only wins in a specific medium-size region and for
           huge products when using many threads. */
        ((len2 >= 8 && len2 <= 75 && bits1 + bits2 >= 800 && bits1 + bits2 <= 4000) ||
            (len1 + len2 >= 5000 && bits1 + bits2 >= 5000 + (len1 + len2) / 10 && flint_get_num_threads() >= 4))
        _fmpz_poly_mul_SS(res, poly1, len1, poly2, len2);
    else
        _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);

#else

    if (len2 < 7)
    {
        _fmpz_poly_mul_classical(res, poly1, len1, poly2, len2);
    }
    else
    {
        mp_size_t limbs1, limbs2;

        limbs1 = (bits1 + FLINT_BITS - 1) / FLINT_BITS;
        limbs2 = (bits2 + FLINT_BITS - 1) / FLINT_BITS;

        if (len1 < 16 && (limbs1 > 12 || limbs2 > 12))
            _fmpz_poly_mul_karatsuba(res, poly1, len1, poly2, len2);
        else if (limbs1 + limbs2 <= 8)
            _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
        else if ((limbs1+limbs2)/2048 > len1 + len2)
            _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
        else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
           _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
        else
           _fmpz_poly_mul_SS(res, poly1, len1, poly2, len2);
    }
#endif
}

void
fmpz_poly_mul(fmpz_poly_t res,
              const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong rlen;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    rlen = len1 + len2 - 1;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        if (len1 >= len2)
            _fmpz_poly_mul(t->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2);
        else
            _fmpz_poly_mul(t->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        if (len1 >= len2)
            _fmpz_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2);
        else
            _fmpz_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1);
    }

    _fmpz_poly_set_length(res, rlen);
}
