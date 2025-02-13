/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_poly.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz * res, const fmpz_t e, const fmpz * f,
                                    slong lenf, const fmpz* finv, slong lenfinv,
                                    const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_powmod_x_fmpz_preinv(res, e, f, lenf, finv, lenfinv, gr_ctx));
}


void
fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res, const fmpz_t e,
                           const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                      const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_poly_powmod_x_fmpz_preinv((gr_poly_struct *) res,
        e,
        (const gr_poly_struct *) f,
        (const gr_poly_struct *) finv, gr_ctx));
}
