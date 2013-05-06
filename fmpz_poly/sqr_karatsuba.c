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
#include "ulong_extras.h"

/*
    For documentation, see fmpz_poly/mul_karatsuba.c
 */

extern void revbin1(fmpz * out, const fmpz * in, long len, long bits);

extern void revbin2(fmpz * out, const fmpz * in, long len, long bits);

extern void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, long bits);

void _fmpz_poly_sqr_kara_recursive(fmpz * out, fmpz * rev,
                                   fmpz * temp, long bits)
{
    long len = (1L << bits);
    long m = len / 2;

    if (len == 1)
    {
        fmpz_mul(out, rev, rev);
        fmpz_zero(out + 1);
        return;
    }

    _fmpz_vec_add(temp, rev, rev + m, m);

    _fmpz_poly_sqr_kara_recursive(out, rev, temp + 2 * m, bits - 1);

    _fmpz_poly_sqr_kara_recursive(out + len, temp, temp + 2 * m, bits - 1);

    _fmpz_poly_sqr_kara_recursive(temp, rev + m, temp + 2 * m, bits - 1);

    _fmpz_vec_sub(out + len, out + len, out, len);
    _fmpz_vec_sub(out + len, out + len, temp, len);

    _fmpz_vec_add_rev(out, temp, bits);
}

void _fmpz_poly_sqr_karatsuba(fmpz * res, const fmpz * poly, long len)
{
    fmpz *rev, *out, *temp;
    long length, loglen = 0;

    if (len == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    while ((1L << loglen) < len)
        loglen++;
    length = (1L << loglen);

    rev  = flint_calloc(3 * length, sizeof(fmpz *));
    out  = rev + length;
    temp = _fmpz_vec_init(2 * length);

    revbin1(rev, poly, len, loglen);

    _fmpz_poly_sqr_kara_recursive(out, rev, temp, loglen);

    _fmpz_vec_zero(res, 2 * len - 1);
    revbin2(res, out, 2 * len - 1, loglen + 1);

    _fmpz_vec_clear(temp, 2 * length);
    flint_free(rev);
}

void fmpz_poly_sqr_karatsuba(fmpz_poly_t res, const fmpz_poly_t poly)
{
    long len;

    if (poly->length == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    len = 2 * poly->length - 1;

    fmpz_poly_fit_length(res, len);

    _fmpz_poly_sqr_karatsuba(res->coeffs, poly->coeffs, poly->length);

    _fmpz_poly_set_length(res, len);
}

