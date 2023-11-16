/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_fmpz_mpoly_q_evaluate(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
{
    ca_t t, u;
    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_fmpz_mpoly_evaluate(t, fmpz_mpoly_q_numref(f), x, mctx, ctx);
    ca_fmpz_mpoly_evaluate(u, fmpz_mpoly_q_denref(f), x, mctx, ctx);
    ca_div(res, t, u, ctx);
    ca_clear(t, ctx);
    ca_clear(u, ctx);
}

void
ca_fmpz_mpoly_q_evaluate_no_division_by_zero(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx)
{
    ca_t t, u;
    ca_init(t, ctx);
    ca_init(u, ctx);
    ca_fmpz_mpoly_evaluate(t, fmpz_mpoly_q_numref(f), x, mctx, ctx);
    ca_fmpz_mpoly_evaluate(u, fmpz_mpoly_q_denref(f), x, mctx, ctx);
    /* todo: write a div function for this */
    ca_inv_no_division_by_zero(u, u, ctx);
    ca_mul(res, t, u, ctx);
    ca_clear(t, ctx);
    ca_clear(u, ctx);
}
