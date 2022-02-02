/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_sqrlow_tiny1(fmpz * res, const fmpz * poly, slong len, slong n)
{
    slong i, j, c;

    _fmpz_vec_zero(res, n);

    for (i = 0; i < len; i++)
    {
        c = poly[i];

        if (c != 0)
        {
            if (2 * i < n)
                res[2 * i] += c * c;

            c *= 2;

            for (j = i + 1; j < FLINT_MIN(len, n - i); j++)
                res[i + j] += poly[j] * c;
        }
    }
}

void _fmpz_poly_sqrlow_tiny2(fmpz * res, const fmpz * poly, slong len, slong n)
{
    slong i, j, k, c, d;
    mp_limb_t hi, lo;
    mp_ptr tmp;
    TMP_INIT;

    TMP_START;

    tmp = TMP_ALLOC(2 * n * sizeof(mp_limb_t));

    flint_mpn_zero(tmp, 2 * n);

    for (i = 0; i < len; i++)
    {
        c = poly[i];

        if (c != 0)
        {
            if (2 * i < n)
            {
                smul_ppmm(hi, lo, c, c);
                add_ssaaaa(tmp[4 * i + 1], tmp[4 * i],
                           tmp[4 * i + 1], tmp[4 * i], hi, lo);
            }

            c *= 2;  /* does not overflow */

            for (j = i + 1; j < FLINT_MIN(len, n - i); j++)
            {
                k = i + j;

                d = poly[j];

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

void _fmpz_poly_sqrlow(fmpz * res, const fmpz * poly, slong len, slong n)
{
    mp_size_t limbs;
    slong bits, rbits;

    len = FLINT_MIN(len, n);

    if (len == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    bits = FLINT_ABS(bits);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX &&
        (len < 50 + 2 * bits || (4 * len >= 3 * n && n < 140 + 6 * bits)))
    {
        rbits = 2 * bits + FLINT_BIT_COUNT(len);

        if (rbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_poly_sqrlow_tiny1(res, poly, len, n);
            return;
        }
        else if (rbits <= 2 * FLINT_BITS - 1)
        {
            _fmpz_poly_sqrlow_tiny2(res, poly, len, n);
            return;
        }
    }

    if (n < 7)
    {
        _fmpz_poly_sqrlow_classical(res, poly, len, n);
        return;
    }

    limbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

    if (n < 16 && limbs > 12)
    {
        int i;
        fmpz *copy;

        copy = flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < len; i++)
            copy[i] = poly[i];
        flint_mpn_zero((mp_ptr) copy + len, n - len);

        _fmpz_poly_sqrlow_karatsuba_n(res, copy, n);

        flint_free(copy);
    }
    else if (limbs <= 4)
        _fmpz_poly_sqrlow_KS(res, poly, len, n);
    else if (limbs/2048 > len)
        _fmpz_poly_sqrlow_KS(res, poly, len, n);
    else if (limbs*FLINT_BITS*4 < len)
       _fmpz_poly_sqrlow_KS(res, poly, len, n);
    else
       _fmpz_poly_mullow_SS(res, poly, len, poly, len, n);
}

void fmpz_poly_sqrlow(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    const slong len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_sqrlow(t, poly, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    n = FLINT_MIN(2 * len - 1, n);

    fmpz_poly_fit_length(res, n);
    _fmpz_poly_sqrlow(res->coeffs, poly->coeffs, len, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
