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

void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz *poly, slong len,
                                  const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);

    if (fmpz_sgn(a) >= 0 && fmpz_cmp(a, fmpz_mod_ctx_modulus(ctx)) < 0)
    {
        GR_MUST_SUCCEED(_gr_poly_evaluate_horner(res, poly, len, a, gr_ctx));
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mod_set_fmpz(t, a, ctx);
        GR_MUST_SUCCEED(_gr_poly_evaluate_horner(res, poly, len, t, gr_ctx));
        fmpz_clear(t);
    }
}

void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz_mod_poly_t poly,
                                      const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    if (res == a)
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_mod_poly_evaluate_fmpz(t, poly->coeffs, poly->length, a, ctx);
        fmpz_swap(res, t);
        fmpz_clear(t);
    }
    else
    {
        _fmpz_mod_poly_evaluate_fmpz(res, poly->coeffs, poly->length, a, ctx);
    }
}
