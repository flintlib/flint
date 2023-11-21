/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_binomial_transform_basecase(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec)
{
    slong n, k;

    fmpz_t t;
    fmpz_init(t);

    for (n = 0; n < len; n++)
    {
        acb_zero(b + n);

        for (k = 0; k < FLINT_MIN(n + 1, alen); k++)
        {
            if (k == 0)
            {
                fmpz_one(t);
            }
            else
            {
                fmpz_mul_si(t, t, -(n - k + 1));
                fmpz_divexact_ui(t, t, k);
            }

            acb_addmul_fmpz(b + n, a + k, t, prec);
        }
    }

    fmpz_clear(t);
}

void
acb_poly_binomial_transform_basecase(acb_poly_t b, const acb_poly_t a, slong len, slong prec)
{
    if (len == 0 || a->length == 0)
    {
        acb_poly_zero(b);
        return;
    }

    if (b == a)
    {
        acb_poly_t c;
        acb_poly_init2(c, len);
        _acb_poly_binomial_transform_basecase(c->coeffs, a->coeffs, a->length, len, prec);
        acb_poly_swap(b, c);
        acb_poly_clear(c);
    }
    else
    {
        acb_poly_fit_length(b, len);
        _acb_poly_binomial_transform_basecase(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _acb_poly_set_length(b, len);
    _acb_poly_normalise(b);
}

