/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

slong fmpz_poly_num_real_roots_sturm(fmpz_poly_t pol)
{
    fmpz_poly_t p0, p1;
    slong i;
    ulong d;
    fmpz_t c;
    int s, s0a, s0b, t;

    if (fmpz_poly_is_zero(pol))
    {
        fprintf(stderr, "ERROR (fmpz_poly_num_real_roots_sturm): zero polynomial\n");
        flint_abort();
    }

    fmpz_init(c);
    fmpz_poly_init(p0);
    fmpz_poly_init(p1);

    fmpz_poly_set(p0, pol);
    fmpz_poly_derivative(p1, p0);

    s0a = s0b = fmpz_sgn(pol->coeffs + pol->length - 1);
    if (pol->length % 2 == 0)
        s0b = -s0b;

    t = 0;
    while (!fmpz_poly_is_zero(p1))
    {
        /* sign change at + infinity */
        s = fmpz_sgn(p1->coeffs + p1->length - 1);
        if (s != s0a)
        {
            t--;
            s0a = s;
        }

        /* sign change at - infinity */
        if (p1->length % 2 == 0) s = -s;
        if (s != s0b)
        {
            t++;
            s0b = s;
        }

        fmpz_poly_swap(p0, p1);
        fmpz_poly_pseudo_rem(p1, &d, p1, p0);
        if ((d%2 == 0) || (fmpz_sgn(p0->coeffs + p0->length - 1) == 1))
            fmpz_poly_neg(p1, p1);
        fmpz_poly_content(c, p1);
        if (!fmpz_is_one(c))
        {
            for (i = 0; i < p1->length; i++)
                fmpz_divexact(p1->coeffs + i, p1->coeffs + i, c);
        }
    }

    fmpz_poly_clear(p0);
    fmpz_poly_clear(p1);
    fmpz_clear(c);

    return t;
}
