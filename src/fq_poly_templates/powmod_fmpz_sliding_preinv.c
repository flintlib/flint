/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e, ulong k,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_powmod_fmpz_sliding_preinv(res, poly, e, k, f, lenf, finv, lenfinv, gr_ctx));
}


void
TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (TEMPLATE(T, poly_t) res,
                                              const TEMPLATE(T, poly_t) poly,
                                              const fmpz_t e, ulong k,
                                              const TEMPLATE(T, poly_t) f,
                                              const TEMPLATE(T, poly_t) finv,
                                              const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_poly_powmod_fmpz_sliding_preinv((gr_poly_struct *) res,
            (const gr_poly_struct *) poly, e, k,
            (const gr_poly_struct *) f,
            (const gr_poly_struct *) finv, gr_ctx));
}


#endif
