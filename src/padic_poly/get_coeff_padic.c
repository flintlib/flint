/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"
#include "padic_poly.h"

void padic_poly_get_coeff_padic(padic_t c, const padic_poly_t poly, slong n,
                                const padic_ctx_t ctx)
{
    if (n < poly->length && !fmpz_is_zero(poly->coeffs + n))
    {
        fmpz_set(padic_unit(c), poly->coeffs + n);
        padic_val(c) = poly->val;
        padic_prec(c) = poly->N;
        padic_reduce(c, ctx);
    }
    else
    {
        padic_zero(c);
    }
}
