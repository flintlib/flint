/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_poly.h"
#include "acb_poly.h"

int
acb_poly_contains_fmpq_poly(const acb_poly_t poly1, const fmpq_poly_t poly2)
{
    slong i;
    fmpq_t t;

    if (poly2->length > poly1->length)
        return 0;

    fmpq_init(t);

    for (i = 0; i < poly2->length; i++)
    {
        fmpq_poly_get_coeff_fmpq(t, poly2, i);
        if (!acb_contains_fmpq(poly1->coeffs + i, t))
        {
            fmpq_clear(t);
            return 0;
        }
    }

    fmpq_clear(t);

    for (i = poly2->length; i < poly1->length; i++)
        if (!acb_contains_zero(poly1->coeffs + i))
            return 0;

    return 1;
}
