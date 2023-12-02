/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int _fmpz_poly_divides(fmpz * q, const fmpz * a,
                       slong len1, const fmpz * b, slong len2)
{
    fmpz * r;

    FLINT_ASSERT(len1 >= len2);
    FLINT_ASSERT(len2 > 0);

    if (!fmpz_divisible(a + 0, b + 0))
        return 0;

    /* heuristic test: see if polys evaluated at 1 divide */
    if (len1 > 1)
    {
        slong i;
        fmpz_t asum, bsum;
        int divisible = 0;

	fmpz_init(asum);
	fmpz_init(bsum);

	for (i = 0; i < len1; i++)
           fmpz_add(asum, asum, a + i);

	for (i = 0; i < len2; i++)
           fmpz_add(bsum, bsum, b + i);

	divisible = fmpz_divisible(asum, bsum);

	fmpz_clear(asum);
	fmpz_clear(bsum);

	if (!divisible)
           return 0;
    }

    r = _fmpz_vec_init(len1);

    if (!_fmpz_poly_divrem(q, r, a, len1, b, len2, 1))
    {
        _fmpz_vec_clear(r, len1);
        return 0;
    }

    FMPZ_VEC_NORM(r, len1);

    _fmpz_vec_clear(r, len1);

    return (len1 == 0);
}

int fmpz_poly_divides(fmpz_poly_t q, const fmpz_poly_t a, const fmpz_poly_t b)
{
    if (fmpz_poly_is_zero(b))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_divides). Division by zero.\n");
    }
    if (fmpz_poly_is_zero(a))
    {
        fmpz_poly_zero(q);
        return 1;
    }
    if (a->length < b->length)
    {
        return 0;
    }

    {
        const slong lenQ = a->length - b->length + 1;
        int res;

        if (q == a || q == b)
        {
            fmpz_poly_t t;

            fmpz_poly_init2(t, lenQ);
            res = _fmpz_poly_divides(t->coeffs, a->coeffs, a->length, b->coeffs, b->length);
            _fmpz_poly_set_length(t, lenQ);
            _fmpz_poly_normalise(t);
            fmpz_poly_swap(q, t);
            fmpz_poly_clear(t);
        }
        else
        {
            fmpz_poly_fit_length(q, lenQ);
            res = _fmpz_poly_divides(q->coeffs, a->coeffs, a->length, b->coeffs, b->length);
            _fmpz_poly_set_length(q, lenQ);
            _fmpz_poly_normalise(q);
        }
        return res;
    }
}
