/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

static void fmpz_poly_evaluate_at_one(fmpz_t res, const fmpz * p, slong len)
{
    _fmpz_vec_sum(res, p, len);
}

slong fmpz_poly_num_real_roots_0_1_sturm(const fmpz_poly_t pol)
{
    fmpz_poly_t p0, p1;
    ulong d;
    fmpz_t c;
    int s, s0a, s0b, t;
    int zero_at_0, zero_at_1;

    if (fmpz_poly_is_zero(pol))
        flint_throw(FLINT_ERROR, "(fmpz_poly_num_real_roots_sturm): zero polynomial\n");

    fmpz_init(c);
    fmpz_poly_init(p0);
    fmpz_poly_init(p1);

    fmpz_poly_set(p0, pol);
    fmpz_poly_derivative(p1, p0);

    s0a = fmpz_sgn(pol->coeffs);
    fmpz_poly_evaluate_at_one(c, pol->coeffs, pol->length);
    s0b = fmpz_sgn(c);

    zero_at_0 = (s0a == 0);
    zero_at_1 = (s0b == 0);

    /* Sturm's theorem counts roots in (0, 1], so count if there's one at 0 */
    t = zero_at_0;

    while (!fmpz_poly_is_zero(p1))
    {
        /* sign change at 0 */
        s = fmpz_sgn(p1->coeffs);
        if (s != 0)
        {
            if (s0a != 0 && s != s0a)
                t++;
            s0a = s;
        }

        /* sign change at 1 */
        fmpz_poly_evaluate_at_one(c, p1->coeffs, p1->length);
        s = fmpz_sgn(c);
        if (s != 0)
        {
            if (s0b != 0 && s != s0b)
                t--;
            s0b = s;
        }

        fmpz_poly_swap(p0, p1);
        fmpz_poly_pseudo_rem(p1, &d, p1, p0);
        if ((d%2 == 0) || (fmpz_sgn(p0->coeffs + p0->length - 1) == 1))
            fmpz_poly_neg(p1, p1);
        fmpz_poly_content(c, p1);
        if (!fmpz_is_one(c))
            _fmpz_vec_scalar_divexact_fmpz(p1->coeffs, p1->coeffs, p1->length, c);
    }

    fmpz_poly_clear(p0);
    fmpz_poly_clear(p1);
    fmpz_clear(c);

    /* We return the count on (0, 1) for consistency with VCA */
    t -= zero_at_1;

    return t;
}

