/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

void _fmpz_mod_poly_pow(fmpz *res, const fmpz *poly, slong len, ulong e,
            const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_pow_ui_binexp(res, poly, len, e, gr_ctx));
}

void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, ulong e,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len = op->length;
    slong rlen;

    if ((len < 2) || (e < UWORD(3)))
    {
        if (e == UWORD(0))
            fmpz_mod_poly_set_ui(rop, 1, ctx);
        else if (len == 0)
            fmpz_mod_poly_zero(rop, ctx);
        else if (len == 1)
        {
            fmpz_mod_poly_fit_length(rop, 1, ctx);
            fmpz_mod_pow_ui(rop->coeffs, op->coeffs, e, ctx);
            _fmpz_mod_poly_set_length(rop, 1);
            _fmpz_mod_poly_normalise(rop);
        }
        else if (e == UWORD(1))
            fmpz_mod_poly_set(rop, op, ctx);
        else  /* e == UWORD(2) */
            fmpz_mod_poly_sqr(rop, op, ctx);
        return;
    }

    rlen = (slong) e * (len - 1) + 1;

    if (rop != op)
    {
        fmpz_mod_poly_fit_length(rop, rlen, ctx);
        _fmpz_mod_poly_pow(rop->coeffs, op->coeffs, len, e, ctx);
        _fmpz_mod_poly_set_length(rop, rlen);
    }
    else
    {
        fmpz *t = _fmpz_vec_init(rlen);

        _fmpz_mod_poly_pow(t, op->coeffs, len, e, ctx);

        _fmpz_vec_clear(rop->coeffs, rop->alloc);
        rop->coeffs = t;
        rop->alloc  = rlen;
        rop->length = rlen;
    }

    _fmpz_mod_poly_normalise(rop);
}
