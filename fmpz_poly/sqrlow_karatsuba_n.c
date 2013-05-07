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

    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    For documentation, see fmpz_poly/mullow_karatsuba_n.c
 */

void _fmpz_poly_sqrlow_kara_recursive(fmpz * out, 
                                      const fmpz * pol, fmpz * temp, long len)
{
    long m1 = len / 2;
    long m2 = len - m1;
    int odd = (len & 1);

    if (len <= 6)
    {
        _fmpz_poly_sqrlow_classical(out, pol, len, len);
        return;
    }

    _fmpz_vec_add(temp + m2, pol, pol + m1, m1);
    if (odd)
        fmpz_set(temp + m2 + m1, pol + 2 * m1);

    _fmpz_poly_sqrlow_kara_recursive(temp, temp + m2, temp + 2 * m2, m2);
    _fmpz_poly_sqrlow_kara_recursive(temp + m2, pol + m1, temp + 2 * m2, m2);

    _fmpz_poly_sqr_karatsuba(out, pol, m1);
    fmpz_zero(out + 2 * m1 - 1);

    _fmpz_vec_sub(temp, temp, out, m2);
    _fmpz_vec_sub(temp, temp, temp + m2, m2);

    if (odd)
        fmpz_set(out + 2 * m1, temp + m2);
    _fmpz_vec_add(out + m1, out + m1, temp, m2);
}

/* 
    Assumes poly1 and poly2 are not length 0.

    We almost get away with temporary space of length 2 * len, 
    but in the recursion we might need 4 * \ceil{len/2}, which 
    exceeds 2 * len by at most 2.
*/
void _fmpz_poly_sqrlow_karatsuba_n(fmpz * res, const fmpz * poly, long n)
{
    fmpz *temp;
    long len, loglen = 0;

    if (n == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    while ((1L << loglen) < n)
        loglen++;
    len = (1L << loglen);

    temp = _fmpz_vec_init(2 * len + 2);

    _fmpz_poly_sqrlow_kara_recursive(res, poly, temp, n);

    _fmpz_vec_clear(temp, 2 * len + 2);
}

void
fmpz_poly_sqrlow_karatsuba_n(fmpz_poly_t res, const fmpz_poly_t poly, long n)
{
    const long len = FLINT_MIN(poly->length, n);
    long i, lenr;

    int clear = 0;
    fmpz *copy;

    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    lenr = 2 * len - 1;
    if (n > lenr)
        n = lenr;

    if (len >= n)
        copy = poly->coeffs;
    else
    {
        copy = flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < len; i++)
            copy[i] = poly->coeffs[i];
        flint_mpn_zero((mp_ptr) copy + len, n - len);
        clear = 1;
    }

    if (res != poly)
    {
        fmpz_poly_fit_length(res, n);
        _fmpz_poly_sqrlow_karatsuba_n(res->coeffs, copy, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_sqrlow_karatsuba_n(t->coeffs, copy, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);

    if (clear)
        flint_free(copy);
}

