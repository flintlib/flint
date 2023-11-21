/*
    Copyright (C) 2012, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "gr_special.h"

void arb_fib_fmpz(arb_t f, const fmpz_t n, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(gr_generic_fib_fmpz(f, n, ctx));
}

void arb_fib_ui(arb_t f, ulong n, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(gr_generic_fib_ui(f, n, ctx));
}
