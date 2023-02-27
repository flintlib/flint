/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void padic_poly_get_coeff_padic(padic_t x, const padic_poly_t f, slong n, 
                                const padic_ctx_t ctx)
{
    if (n < f->length && !fmpz_is_zero(f->coeffs + n))
    {
        fmpz_set(padic_unit(x), f->coeffs + n);
        padic_val(x) = f->val;
        padic_prec(x) = f->N;
        padic_reduce(x, ctx);
    }
    else
    {
        padic_zero(x);
    }
}

