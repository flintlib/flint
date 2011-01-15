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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    The following is a port of the FLINT1 implementation, although the 
    handling of the base cases has been changed.
 */
static void 
_yyy(fmpz_t output, const fmpz * poly, long len, const fmpz_t val)
{
    long i, log, lenT;
    fmpz *half, *temp;
    fmpz *coeff_h, *coeff_t;
    fmpz_t pow;

    if (len < 3)
    {
        _fmpz_poly_evaluate_horner_fmpz(output, poly, len, val);
        return;
    }

    fmpz_init(pow);

    temp = _fmpz_vec_init((len + 1) / 2);

    coeff_t = (fmpz *) poly;
    coeff_h = temp;

    for (i = 0; i < len / 2; i++)
    {
        fmpz_mul(coeff_h, coeff_t + 1, val);
        fmpz_add(coeff_h, coeff_h, coeff_t);
        coeff_t += 2;
        coeff_h ++;
    }

    if (len & 1L)
        fmpz_set(coeff_h, coeff_t);

    fmpz_mul(pow, val, val);

    for (lenT = (len + 1) / 2, log = 1; lenT > 2; lenT = (lenT + 1) / 2, log++)
    {
        half = _fmpz_vec_init((lenT + 1) / 2);

        coeff_t = temp;
        coeff_h = half;

        for (i = 0; i < lenT / 2; i++)
        {
            fmpz_mul(coeff_h, coeff_t + 1, pow);
            fmpz_add(coeff_h, coeff_h, coeff_t);
            coeff_t += 2;
            coeff_h ++;
        }

        fmpz_mul(pow, pow, pow);

        if (lenT & 1L)
            fmpz_set(coeff_h, coeff_t);

        _fmpz_vec_swap(temp, half, FLINT_MIN(lenT, (lenT + 1)/2));
        _fmpz_vec_clear(half, (lenT + 1) / 2);
    }

    fmpz_mul(output, temp + 1, pow);
    fmpz_add(output, output, temp);

    _fmpz_vec_clear(temp, (len + 1) / 2);
    fmpz_clear(pow);
}

/*
    Computes res = (poly, len)(x) where y is an array such that 
    y[i] = x^{2^i}.  Allows zero-padding.
 */
static void 
_xxx_recursive(fmpz_t res, const fmpz * poly, long len, const fmpz_t x, 
                           const fmpz * y)
{
    if (len == 1)
    {
        fmpz_set(res, poly);
    }
    else
    {
        long k;
        fmpz *x0 = res;
        fmpz_t x1;

        k = FLINT_BIT_COUNT(len - 1);

        fmpz_init(x1);
        _xxx_recursive(x0, poly,                   1L << (k - 1),         x, y);
        _xxx_recursive(x1, poly + (1L << (k - 1)), len - (1L << (k - 1)), x, y);
        fmpz_addmul(x0, y + (k - 1), x1);
        fmpz_clear(x1);
    }
}

static void 
_xxx(fmpz_t res, const fmpz * poly, long len, const fmpz_t x)
{
    if (len == 1)
    {
        fmpz_set(res, poly);
    }
    else
    {
        long i, k;
        fmpz *y;

        k = FLINT_BIT_COUNT(len - 1);  /* 2^{k-1} < len <= 2^k */
        y = _fmpz_vec_init(2 * k);     /* x^{2^0}, x^{2^1}, ..., x^{2^{k-1}} */

        *y = *x;
        for (i = 1; i < k; i++)
            fmpz_mul(y + i, y + (i - 1), y + (i - 1));

        _xxx_recursive(res, poly, len, x, y);

        *y = 0L;
        _fmpz_vec_clear(y, k);
    }
}

void 
_fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz * poly, long len, 
                                                const fmpz_t x)
{
    long c, h, i, k;
    fmpz *y, *T, *t = res, *u;

    h = FLINT_BIT_COUNT(len - 1);  /* 2^{h-1} < len <= 2^h */
    y = _fmpz_vec_init(2 * h + 2); /* x^{2^0}, x^{2^1}, ..., x^{2^{h-1}} */
    T = y + h;
    u = y + 2 * h + 1;

    *y = *x;
    for (i = 1; i < h; i++)
        fmpz_mul(y + i, y + (i - 1), y + (i - 1));

    for (i = 0; i < len; )
    {
        fmpz_set(t, poly + i);
        ++i;
        count_trailing_zeros(c, i);
        for (k = 0; k < c; k++)
        {
            fmpz_mul(u, y + k, t);
            fmpz_add(t, T + k, u);
        }
        fmpz_swap(T + k, t);
    }
    fmpz_swap(t, T + k);

    for ( ; k < h; k++)
    {
        if ((len - 1) & (1L << k))
        {
            fmpz_mul(u, y + k, t);
            fmpz_add(t, T + k, u);
        }
    }

    *y = 0L;
    _fmpz_vec_clear(y, 2 * h + 2);
}

void
fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz_poly_t poly, 
                                   const fmpz_t a)
{
    if (fmpz_poly_is_zero(poly))
    {
        fmpz_zero(res);
        return;
    }

    if (res == a)
    {
        fmpz_t t;

        fmpz_init(t);
        _fmpz_poly_evaluate_divconquer_fmpz(t, poly->coeffs, poly->length, a);
        fmpz_swap(res, t);
        fmpz_clear(t);
    }
    else
        _fmpz_poly_evaluate_divconquer_fmpz(res, poly->coeffs, poly->length, a);
}

