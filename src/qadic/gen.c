/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "qadic.h"

void qadic_gen(qadic_t x, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(x);
    const slong d = qadic_ctx_degree(ctx);

    if (d > 1)
    {
        if (N > 0)
        {
            padic_poly_fit_length(x, 2);
            fmpz_zero(x->coeffs + 0);
            fmpz_one(x->coeffs + 1);
            _padic_poly_set_length(x, 2);
            x->val = 0;
        }
        else
        {
            padic_poly_zero(x);
        }
    }
    else
    {
        padic_poly_fit_length(x, 1);
	fmpz_neg(x->coeffs + 0, ctx->a + 0);
	_padic_poly_set_length(x, 1);
	x->val = 0;
    }
}
