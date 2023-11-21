/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "gr_special.h"

void
arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)
{
    gr_ctx_t ctx;
    fmpz_t m;
    gr_ctx_init_real_arb(ctx, prec);
    fmpz_init_set_ui(m, n);
    GR_MUST_SUCCEED(gr_generic_chebyshev_u2_fmpz(a, b, m, x, ctx));
    fmpz_clear(m);
}

void
arb_chebyshev_u_ui(arb_t y, ulong n, const arb_t x, slong prec)
{
    gr_ctx_t ctx;
    fmpz_t m;
    gr_ctx_init_real_arb(ctx, prec);
    fmpz_init_set_ui(m, n);
    GR_MUST_SUCCEED(gr_generic_chebyshev_u_fmpz(y, m, x, ctx));
    fmpz_clear(m);
}
