/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

void _fmpz_mod_poly_derivative(fmpz *res, const fmpz *poly, slong len,
                               const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_derivative(res, poly, len, gr_ctx));
}

void fmpz_mod_poly_derivative(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len = poly->length;

    if (len < 2)
    {
        fmpz_mod_poly_zero(res, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, len - 1, ctx);
        _fmpz_mod_poly_derivative(res->coeffs, poly->coeffs, len, ctx);
        _fmpz_mod_poly_set_length(res, len - 1);
        _fmpz_mod_poly_normalise(res);
    }
}
