/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

void padic_poly_neg(padic_poly_t f, const padic_poly_t g, 
                    const padic_ctx_t ctx)
{
    const slong len = g->length;

    if (len == 0 || g->val >= padic_poly_prec(f))
    {
        padic_poly_zero(f);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        padic_poly_fit_length(f, len);
        _padic_poly_set_length(f, len);
        f->val = g->val;

        alloc = _padic_ctx_pow_ui(pow, padic_poly_prec(f) - f->val, ctx);

        if (padic_poly_prec(f) >= padic_poly_prec(g))  /* No reduction */
        {
            _fmpz_mod_poly_neg(f->coeffs, g->coeffs, len, pow);
        }
        else  /* Reduction necessary? */
        {
            _fmpz_vec_scalar_mod_fmpz(f->coeffs, g->coeffs, len, pow);
            _fmpz_mod_poly_neg(f->coeffs, f->coeffs, len, pow);
            _padic_poly_normalise(f);
        }

        if (alloc)
            fmpz_clear(pow);
    }
}

