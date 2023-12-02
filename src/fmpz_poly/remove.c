/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

slong
fmpz_poly_remove(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    fmpz_poly_t p, q;
    fmpz_t p1sum, p2sum, qsum;
    slong i;

    if (poly2->length == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_remove): Division by zero.\n");
    }

    if (poly2->length == 1 && fmpz_is_pm1(poly2->coeffs + 0))
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_remove): Divisor must not be a unit.\n");
    }

    if (poly2->length > poly1->length)
    {
        fmpz_poly_set(res, poly1);
        return 0;
    }

    fmpz_init(p1sum);
    fmpz_init(p2sum);
    fmpz_init(qsum);

    for (i = 0; i < poly1->length; i++)
        fmpz_add(p1sum, p1sum, poly1->coeffs + i);

    for (i = 0; i < poly2->length; i++)
        fmpz_add(p2sum, p2sum, poly2->coeffs + i);

    fmpz_abs(p1sum, p1sum);
    fmpz_abs(p2sum, p2sum);

    if (fmpz_is_zero(p2sum))
    {
        if (!fmpz_is_zero(p1sum))
        {
            fmpz_poly_set(res, poly1);
            i = 0;
            goto cleanup;
        } else
            i = (poly1->length - 1)/(poly2->length - 1);
    } else if (fmpz_is_zero(p1sum) || fmpz_is_one(p2sum))
        i = (poly1->length - 1)/(poly2->length - 1);
    else
        i = fmpz_remove(qsum, p1sum, p2sum);

    if (i > 0)
    {
        fmpz_poly_init(q);
        fmpz_poly_init(p);

        fmpz_poly_pow(p, poly2, i);

        while (i > 0 && !fmpz_poly_divides(q, poly1, p))
        {
            fmpz_poly_div(p, p, poly2);
            i--;
        }

        if (i == 0)
            fmpz_poly_set(res, poly1);
        else
            fmpz_poly_set(res, q);

        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    } else
        fmpz_poly_set(res, poly1);

cleanup:

    fmpz_clear(qsum);
    fmpz_clear(p1sum);
    fmpz_clear(p2sum);

    return i;
}
