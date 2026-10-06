/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "gr.h"
#include "gr_poly.h"

int
_fmpz_mod_poly_is_squarefree(const fmpz * f, slong len, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_struct t;
    truth_t res;

    if (len <= 2)
        return len != 0;

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(gr_ctx, T_TRUE));

    t.coeffs = (fmpz *) f;
    t.length = len;
    t.alloc = len;

    res = gr_poly_is_squarefree(&t, gr_ctx);

    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "Exception in _fmpz_mod_poly_is_squarefree: "
                                 "unable to decide\n");

    return res == T_TRUE;
}

int fmpz_mod_poly_is_squarefree(const fmpz_mod_poly_t f,
                                                      const fmpz_mod_ctx_t ctx)
{
    return _fmpz_mod_poly_is_squarefree(f->coeffs, f->length, ctx);
}
