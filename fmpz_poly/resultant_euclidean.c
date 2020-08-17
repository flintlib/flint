/*
    Copyright (C) 2010 Sebastian Pancratz

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

void
_fmpz_poly_resultant_euclidean(fmpz_t res, const fmpz * poly1, slong len1, 
                                 const fmpz * poly2, slong len2)
{
    if (len2 == 1)
    {
        fmpz_pow_ui(res, poly2, len1 - 1);
    }
    else
    {
        fmpz_t a, b, g, h, t;
        fmpz *A, *B, *W;
        const slong alloc = len1 + len2;
        slong sgn = 1;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(g);
        fmpz_init(h);
        fmpz_init(t);

        A = W = _fmpz_vec_init(alloc);
        B = W + len1;

        _fmpz_poly_content(a, poly1, len1);
        _fmpz_poly_content(b, poly2, len2);
        _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, a);
        _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, b);

        fmpz_one(g);
        fmpz_one(h);

        fmpz_pow_ui(a, a, len2 - 1);
        fmpz_pow_ui(b, b, len1 - 1);
        fmpz_mul(t, a, b);

        do
        {
            const slong d = len1 - len2;

            if (!(len1 & WORD(1)) & !(len2 & WORD(1)))
                sgn = -sgn;

            _fmpz_poly_pseudo_rem_cohen(A, A, len1, B, len2);

            FMPZ_VEC_NORM(A, len1);

            if (len1 == 0)
            {
                fmpz_zero(res);
                goto cleanup;
            }

            {
                fmpz * T;
                slong len;
                T = A, A = B, B = T;
                len = len1, len1 = len2, len2 = len;
            }

            fmpz_pow_ui(a, h, d);
            fmpz_mul(b, g, a);
            _fmpz_vec_scalar_divexact_fmpz(B, B, len2, b);

            fmpz_pow_ui(g, A + (len1 - 1), d);
            fmpz_mul(b, h, g);
            fmpz_divexact(h, b, a);
            fmpz_set(g, A + (len1 - 1));

        } while (len2 > 1);

        fmpz_pow_ui(g, h, len1 - 1);
        fmpz_pow_ui(b, B + (len2 - 1), len1 - 1);
        fmpz_mul(a, h, b);
        fmpz_divexact(h, a, g);

        fmpz_mul(res, t, h);
        if (sgn < 0)
            fmpz_neg(res, res);

      cleanup:

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(g);
        fmpz_clear(h);
        fmpz_clear(t);

        _fmpz_vec_clear(W, alloc);
    }
}

void
fmpz_poly_resultant_euclidean(fmpz_t res, const fmpz_poly_t poly1, 
                                const fmpz_poly_t poly2)
{
    const slong len1 = poly1->length, len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_zero(res);
        return;
    }

    if (len1 >= len2)
        _fmpz_poly_resultant_euclidean(res, poly1->coeffs, len1, poly2->coeffs, len2);
    else
    {
        _fmpz_poly_resultant_euclidean(res, poly2->coeffs, len2, poly1->coeffs, len1);
        if ((len1 > 1) && (!(len1 & WORD(1)) & !(len2 & WORD(1))))
            fmpz_neg(res, res);
    }
}
