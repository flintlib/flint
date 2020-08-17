/*
    Copyright (C) 2010 William Hart

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

/* Assumes poly1 and poly2 are not length 0. */
void
_fmpz_poly_mulhigh_classical(fmpz * res, const fmpz * poly1,
                             slong len1, const fmpz * poly2, slong len2,
                             slong start)
{
    _fmpz_vec_zero(res, start);

    if (len1 == 1 && len2 == 1) /* Special case if the length of both inputs is 1 */
    {
        if (start == 0)
            fmpz_mul(res, poly1, poly2);
    }
    else                        /* Ordinary case */
    {
        slong i, m, n;

        /* Set res[i] = poly1[i]*poly2[0] */
        if (start < len1)
            _fmpz_vec_scalar_mul_fmpz(res + start, poly1 + start, len1 - start,
                                      poly2);

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        m = FLINT_MAX(len1 - 1, start);
        _fmpz_vec_scalar_mul_fmpz(res + m, poly2 + m - len1 + 1,
                                  len2 - 1 + len1 - m, poly1 + len1 - 1);

        /* out[i+j] += in1[i]*in2[j] */
        m = FLINT_MAX(start, len2 - 1);
        for (i = m - len2 + 1; i < len1 - 1; i++)
        {
            n = FLINT_MAX(i + 1, start);
            _fmpz_vec_scalar_addmul_fmpz(res + n, poly2 + n - i, len2 + i - n,
                                         poly1 + i);
        }
    }
}

void
fmpz_poly_mulhigh_classical(fmpz_poly_t res,
                            const fmpz_poly_t poly1, const fmpz_poly_t poly2,
                            slong start)
{
    slong len_out = poly1->length + poly2->length - 1;

    if (poly1->length == 0 || poly2->length == 0 || start >= len_out)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t temp;
        fmpz_poly_init2(temp, len_out);
        _fmpz_poly_mulhigh_classical(temp->coeffs, poly1->coeffs,
                                     poly1->length, poly2->coeffs,
                                     poly2->length, start);
        fmpz_poly_swap(res, temp);
        fmpz_poly_clear(temp);
    }
    else
    {
        fmpz_poly_fit_length(res, len_out);
        _fmpz_poly_mulhigh_classical(res->coeffs, poly1->coeffs, poly1->length,
                                     poly2->coeffs, poly2->length, start);
    }

    _fmpz_poly_set_length(res, len_out);
}
