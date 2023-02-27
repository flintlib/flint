/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void padic_poly_reduce(padic_poly_t poly, const padic_ctx_t ctx)
{
    const slong N = padic_poly_prec(poly);

    if (poly->length > 0)
    {
        if (poly->val >= N)
        {
            padic_poly_zero(poly);
        }
        else
        {
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, N - poly->val, ctx);

            _fmpz_vec_scalar_mod_fmpz(poly->coeffs, poly->coeffs, poly->length, pow);

            if (alloc)
                fmpz_clear(pow);

            _padic_poly_normalise(poly);

            if (poly->length == 0)
            {
                poly->val = 0;
            }
        }
    }
}

