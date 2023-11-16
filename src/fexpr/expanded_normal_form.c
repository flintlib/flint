/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"
#include "fexpr.h"

int
fexpr_expanded_normal_form(fexpr_t res, const fexpr_t expr, ulong flags)
{
    fexpr_vec_t args;
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_q_t frac;
    int success;

    fexpr_vec_init(args, 0);

    fexpr_arithmetic_nodes(args, expr);
    _fexpr_vec_sort_fast(args->entries, args->length);

    /* todo: when length == 0, use fmpq arithmetic instead */
    fmpz_mpoly_ctx_init(ctx, FLINT_MAX(args->length, 1), ORD_LEX);
    fmpz_mpoly_q_init(frac, ctx);

    success = fexpr_get_fmpz_mpoly_q(frac, expr, args, ctx);

    if (success)
        fexpr_set_fmpz_mpoly_q(res, frac, args, ctx);
    else
        fexpr_set(res, expr);

    fmpz_mpoly_q_clear(frac, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    fexpr_vec_clear(args);

    return success;
}
