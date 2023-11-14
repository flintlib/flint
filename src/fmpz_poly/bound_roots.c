/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

/* quotient of absolute value, rounded up */
static inline void
fmpz_cdiv_abs_q(fmpz_t q, const fmpz_t x, const fmpz_t y)
{
    if (fmpz_sgn(x) == fmpz_sgn(y))
    {
        fmpz_cdiv_q(q, x, y);
    }
    else
    {
        fmpz_fdiv_q(q, x, y);
        fmpz_neg(q, q);
    }
}

void
_fmpz_poly_bound_roots(fmpz_t bound, const fmpz * poly, slong len)
{

    if (len <= 1)
    {
        fmpz_zero(bound);
    }
    else if (len == 2)
    {
        fmpz_cdiv_abs_q(bound, poly + 0, poly + 1);
    }
    else
    {
        slong i;
        fmpz_t t;

        fmpz_init(t);

        fmpz_mul_2exp(t, poly + len - 1, 1);
        fmpz_cdiv_abs_q(bound, poly + 0, t);
        fmpz_root(bound, bound, len - 1);

        for (i = 1; i < len - 1; i++)
        {
            fmpz_cdiv_abs_q(t, poly + len - i - 1, poly + len - 1);
            fmpz_root(t, t, i);
            fmpz_add_ui(t, t, 1);

            if (fmpz_cmp(t, bound) > 0)
                fmpz_swap(t, bound);
        }

        fmpz_mul_2exp(bound, bound, 1);

        fmpz_clear(t);
    }
}

void
fmpz_poly_bound_roots(fmpz_t bound, const fmpz_poly_t poly)
{
    _fmpz_poly_bound_roots(bound, poly->coeffs, poly->length);
}

