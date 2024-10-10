/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_q.h"
#include "gr.h"
#include "gr_generic.h"

static void
fmpz_mpoly_evaluate_acb(acb_t res, const fmpz_mpoly_t pol, acb_srcptr x, slong prec, const fmpz_mpoly_ctx_t ctx)
{
    gr_ctx_t CC;
    gr_ctx_init_complex_acb(CC, prec);

    if (pol->length <= 6 && pol->bits <= FLINT_BITS)
    {
        if (gr_fmpz_mpoly_evaluate_iter(res, pol, x, ctx, CC) != GR_SUCCESS)
            acb_indeterminate(res);
    }
    else
    {
        if (gr_fmpz_mpoly_evaluate_horner(res, pol, x, ctx, CC) != GR_SUCCESS)
            acb_indeterminate(res);
    }
}

void
fmpz_mpoly_q_evaluate_acb(acb_t res, const fmpz_mpoly_q_t f, acb_srcptr x, slong prec, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(f), ctx))
    {
        fmpz_mpoly_evaluate_acb(res, fmpz_mpoly_q_numref(f), x, prec, ctx);
    }
    else
    {
        acb_t t;
        acb_init(t);
        fmpz_mpoly_evaluate_acb(t, fmpz_mpoly_q_denref(f), x, prec, ctx);
        if (acb_contains_zero(t))
        {
            acb_indeterminate(res);
        }
        else
        {
            fmpz_mpoly_evaluate_acb(res, fmpz_mpoly_q_numref(f), x, prec, ctx);
            acb_div(res, res, t, prec);
        }
        acb_clear(t);
    }
}
