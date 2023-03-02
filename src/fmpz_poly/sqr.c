/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_sqr_tiny1(fmpz * res, const fmpz * poly, slong len)
{
    slong i, j, c;

    _fmpz_vec_zero(res, 2 * len - 1);

    for (i = 0; i < len; i++)
    {
        c = poly[i];

        if (c != 0)
        {
            res[2 * i] += c * c;

            c *= 2;

            for (j = i + 1; j < len; j++)
                res[i + j] += poly[j] * c;
        }
    }
}

void _fmpz_poly_sqr_tiny2(fmpz * res, const fmpz * poly, slong len)
{
    slong i, j, k, c, d;
    mp_limb_t hi, lo;
    mp_ptr tmp;
    TMP_INIT;

    TMP_START;

    tmp = TMP_ALLOC(2 * (2 * len - 1) * sizeof(mp_limb_t));

    flint_mpn_zero(tmp, 2 * (2 * len - 1));

    for (i = 0; i < len; i++)
    {
        c = poly[i];

        if (c != 0)
        {
            smul_ppmm(hi, lo, c, c);
            add_ssaaaa(tmp[4 * i + 1], tmp[4 * i],
                       tmp[4 * i + 1], tmp[4 * i], hi, lo);

            c *= 2;  /* does not overflow */

            for (j = i + 1; j < len; j++)
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

    for (i = 0; i < 2 * len - 1; i++)
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

void _fmpz_poly_sqr(fmpz * res, const fmpz * poly, slong len)
{
    mp_size_t limbs;
    slong bits, rbits;

    if (len == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    bits = FLINT_ABS(bits);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX && len < 50 + 3 * bits)
    {
        rbits = 2 * bits + FLINT_BIT_COUNT(len);

        if (rbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_poly_sqr_tiny1(res, poly, len);
            return;
        }
        else if (rbits <= 2 * FLINT_BITS - 1)
        {
            _fmpz_poly_sqr_tiny2(res, poly, len);
            return;
        }
    }

    if (len < 7)
    {
        _fmpz_poly_sqr_classical(res, poly, len);
        return;
    }

    limbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

    if (len < 16 && limbs > 12)
        _fmpz_poly_sqr_karatsuba(res, poly, len);
    else if (limbs <= 4)
        _fmpz_poly_sqr_KS(res, poly, len);
    else if (limbs/2048 > len)
        _fmpz_poly_sqr_KS(res, poly, len);
    else if (limbs*FLINT_BITS*4 < len)
       _fmpz_poly_sqr_KS(res, poly, len);
    else
       _fmpz_poly_mul_SS(res, poly, len, poly, len);
}

void fmpz_poly_sqr(fmpz_poly_t res, const fmpz_poly_t poly)
{
    slong len = poly->length;
    slong rlen;

    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    rlen = 2 * len - 1;

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_sqr(t->coeffs, poly->coeffs, len);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_sqr(res->coeffs, poly->coeffs, len);
    }

    _fmpz_poly_set_length(res, rlen);
}

