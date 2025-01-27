/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

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
_fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, const fmpz * poly3, slong len3,
                              const fmpz * poly3inv, slong len3inv, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_compose_mod_brent_kung_preinv(res, poly1, len1, poly2, poly3, len3, poly3inv, len3inv, gr_ctx));
}

void fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                         const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_poly_compose_mod_brent_kung_preinv((gr_poly_struct *) res,
        (const gr_poly_struct *) poly1,
        (const gr_poly_struct *) poly2,
        (const gr_poly_struct *) poly3,
        (const gr_poly_struct *) poly3inv, gr_ctx));
}
