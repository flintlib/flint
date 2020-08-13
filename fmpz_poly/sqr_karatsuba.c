/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Sebastian Pancratz

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
#include "ulong_extras.h"

/*
    For documentation, see fmpz_poly/mul_karatsuba.c
 */

extern void revbin1(fmpz * out, const fmpz * in, slong len, slong bits);

extern void revbin2(fmpz * out, const fmpz * in, slong len, slong bits);

extern void _fmpz_vec_add_rev(fmpz * in1, fmpz * in2, slong bits);

void _fmpz_poly_sqr_kara_recursive(fmpz * out, fmpz * rev,
                                   fmpz * temp, slong bits)
{
    slong len = (WORD(1) << bits);
    slong m = len / 2;

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

void _fmpz_poly_sqr_karatsuba(fmpz * res, const fmpz * poly, slong len)
{
    fmpz *rev, *out, *temp;
    slong length, loglen = 0;

    if (len == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    while ((WORD(1) << loglen) < len)
        loglen++;
    length = (WORD(1) << loglen);

    rev  = flint_calloc(3 * length, sizeof(fmpz));
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
    slong len;

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

