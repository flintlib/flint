/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "fmpz_mpoly_q.h"

/* stupid algorithm, just to have something working... */
void
fmpz_mpoly_evaluate_acb(acb_t res, const fmpz_mpoly_t pol, acb_srcptr x, slong prec, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, len, nvars;
    acb_t s, t, u;
    ulong * exp; /* todo: support fmpz exponents */

    len = fmpz_mpoly_length(pol, ctx);

    if (len == 0)
    {
        acb_zero(res);
        return;
    }

    if (len == 1 && fmpz_mpoly_is_fmpz(pol, ctx))
    {
        acb_set_round_fmpz(res, pol->coeffs, prec);
        return;
    }

    nvars = ctx->minfo->nvars;
    exp = flint_malloc(sizeof(ulong) * nvars);

    acb_init(s);
    acb_init(t);
    acb_init(u);

    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_get_term_exp_ui(exp, pol, i, ctx);

        acb_one(t);

        for (j = 0; j < nvars; j++)
        {
            if (exp[j] == 1)
            {
                acb_mul(t, t, x + j, prec);
            }
            else if (exp[j] >= 2)
            {
                acb_pow_ui(u, x + j, exp[j], prec);
                acb_mul(t, t, u, prec);
            }
        }

        acb_addmul_fmpz(s, t, pol->coeffs + i, prec);
    }

    acb_swap(res, s);

    flint_free(exp);
    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
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

