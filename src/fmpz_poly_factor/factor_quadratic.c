/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"

void _fmpz_poly_factor_quadratic(fmpz_poly_factor_t fac, 
            const fmpz_poly_t f, slong exp)
{
    fmpz_t D;
    const fmpz *a, *b, *c;

    c = f->coeffs;
    b = f->coeffs + 1;
    a = f->coeffs + 2;

    fmpz_init(D);

    fmpz_mul(D, a, c);
    fmpz_mul_2exp(D, D, 2);
    fmpz_submul(D, b, b);
    fmpz_neg(D, D);

    if (fmpz_is_square(D))
    {
        fmpz_poly_t t;
        fmpz_t g;
        fmpz_poly_init2(t, 2);
        fmpz_init(g);
        _fmpz_poly_set_length(t, 2);

        fmpz_sqrt(D, D);

        fmpz_mul_2exp(t->coeffs + 1, a, 1);
        fmpz_sub(t->coeffs, b, D);
        fmpz_poly_content(g, t);
        fmpz_poly_scalar_divexact_fmpz(t, t, g);

        if (fmpz_is_zero(D))
        {
            fmpz_poly_factor_insert(fac, t, 2 * exp);
        }
        else
        {
            fmpz_poly_factor_insert(fac, t, exp);

            fmpz_mul_2exp(t->coeffs + 1, a, 1);
            fmpz_add(t->coeffs, b, D);
            fmpz_poly_content(g, t);
            fmpz_poly_scalar_divexact_fmpz(t, t, g);
            fmpz_poly_factor_insert(fac, t, exp);
        }

        fmpz_poly_clear(t);
        fmpz_clear(g);
    }
    else
    {
        fmpz_poly_factor_insert(fac, f, exp);
    }

    fmpz_clear(D);
}
