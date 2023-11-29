/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_poly.h"

slong fmpq_poly_remove(fmpq_poly_t q, const fmpq_poly_t poly1,
                                                       const fmpq_poly_t poly2)
{
    fmpq_poly_t p1, p2, qpoly, pow;
    fmpq_t c1, c2;
    fmpz_t p1sum, p2sum, qsum;
    slong i, len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len2 == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_remove): Division by zero.\n");
    }

    if (len2 == 1)
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_remove): Divisor must not be a unit.\n");
    }

    if (len2 > len1)
    {
    	fmpq_poly_set(q, poly1);
	return 0;
    }

    fmpq_poly_init(qpoly);
    fmpq_poly_init(p1);
    fmpq_poly_init(p2);
    fmpq_poly_init(pow);
    fmpq_init(c1);
    fmpq_init(c2);

    fmpq_poly_content(c1, poly1);
    fmpq_poly_content(c2, poly2);

    fmpq_poly_scalar_div_fmpq(p1, poly1, c1);
    fmpq_poly_scalar_div_fmpq(p2, poly2, c2);

    fmpz_init(p1sum);
    fmpz_init(p2sum);
    fmpz_init(qsum);

    for (i = 0; i < poly1->length; i++)
        fmpz_add(p1sum, p1sum, p1->coeffs + i);

    for (i = 0; i < poly2->length; i++)
        fmpz_add(p2sum, p2sum, p2->coeffs + i);

    fmpz_abs(p1sum, p1sum);
    fmpz_abs(p2sum, p2sum);

    if (fmpz_is_zero(p2sum))
    {
        if (!fmpz_is_zero(p1sum))
        {
            fmpq_poly_set(q, poly1);
            i = 0;
            goto cleanup;
        } else
            i = (p1->length - 1)/(p2->length - 1);
    } else if (fmpz_is_zero(p1sum) || fmpz_is_one(p2sum))
        i = (p1->length - 1)/(p2->length - 1);
    else
        i = fmpz_remove(qsum, p1sum, p2sum);

    fmpq_poly_pow(pow, p2, i);

    while (i > 0 && !fmpq_poly_divides(qpoly, p1, pow))
    {
	fmpq_poly_div(pow, pow, p2);
	i--;
    }

    if (i == 0)
       fmpq_poly_set(q, poly1);
    else
    {
        fmpq_pow_si(c2, c2, i);
        fmpq_div(c1, c1, c2);

        fmpq_poly_scalar_mul_fmpq(q, qpoly, c1);
    }

cleanup:

    fmpz_clear(qsum);
    fmpz_clear(p1sum);
    fmpz_clear(p2sum);
    fmpq_clear(c2);
    fmpq_clear(c1);
    fmpq_poly_clear(pow);
    fmpq_poly_clear(p2);
    fmpq_poly_clear(p1);
    fmpq_poly_clear(qpoly);

    return i;
}

