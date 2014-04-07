/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

static __inline__ int
_fmpz_poly_mul_tiny_case(const fmpz * poly1, slong len1,
                         const fmpz * poly2, slong len2)
{
    slong i;
    mp_bitcnt_t bits;
    fmpz c;
    ulong max1, max2;

    max1 = max2 = 0;

    for (i = 0; i < len1; i++)
    {
        c = poly1[i];
        if (COEFF_IS_MPZ(c))
            return 0;
        max1 |= FLINT_ABS(c);
    }

    for (i = 0; i < len2; i++)
    {
        c = poly2[i];
        if (COEFF_IS_MPZ(c))
            return 0;
        max2 |= FLINT_ABS(c);
    }

    bits = FLINT_BIT_COUNT(max1) + FLINT_BIT_COUNT(max2) + FLINT_BIT_COUNT(len2);

    if (bits <= FLINT_BITS - 2)
        return 1;
    else if (bits <= 2 * FLINT_BITS - 1)
        return 2;
    else
        return 0;
}

void
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

void
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
    mp_size_t limbs1, limbs2;

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

    if (len2 < 90)
    {
        int tiny_case = _fmpz_poly_mul_tiny_case(poly1, len1, poly2, len2);

        if (len2 < 50 && tiny_case == 1)
        {
            _fmpz_poly_mul_tiny1(res, poly1, len1, poly2, len2);
            return;
        }

        if (tiny_case == 2)
        {
            _fmpz_poly_mul_tiny2(res, poly1, len1, poly2, len2);
            return;
        }
    }

    if (len2 < 7)
    {
        _fmpz_poly_mul_classical(res, poly1, len1, poly2, len2);
        return;
    }

    limbs1 = _fmpz_vec_max_limbs(poly1, len1);
    limbs2 = _fmpz_vec_max_limbs(poly2, len2);

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
