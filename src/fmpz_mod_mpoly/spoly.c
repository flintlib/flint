/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets
    
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mpoly.h"

void
fmpz_mod_mpoly_spoly(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong n, i;
    ulong * exp, * expf, * expg;
    fmpz_mod_mpoly_t T, U;

    if (f->length == 0 || g->length == 0)
    {
        fmpz_mod_mpoly_zero(res, ctx);
        return;
    }

    n = ctx->minfo->nvars;

    exp = flint_malloc(sizeof(ulong) * n);
    expf = flint_malloc(sizeof(ulong) * n);
    expg = flint_malloc(sizeof(ulong) * n);
    fmpz_mod_mpoly_init(T, ctx);
    fmpz_mod_mpoly_init(U, ctx);

    fmpz_mod_mpoly_get_term_exp_ui(expf, f, 0, ctx);
    fmpz_mod_mpoly_get_term_exp_ui(expg, g, 0, ctx);

    for (i = 0; i < n; i++)
        exp[i] = FLINT_MAX(expf[i], expg[i]);

    for (i = 0; i < n; i++)
    {
        expf[i] = exp[i] - expf[i];
        expg[i] = exp[i] - expg[i];
    }

    if (fmpz_is_one(f->coeffs))
        fmpz_mod_mpoly_set_coeff_ui_ui(T, 1, expf, ctx);
    else
    {
        fmpz_t c;
        fmpz_init(c);
        fmpz_mod_inv(c, f->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_set_coeff_fmpz_ui(T, c, expf, ctx);
        fmpz_clear(c);
    }

    fmpz_mod_mpoly_mul(T, T, f, ctx);
    
    if (fmpz_is_one(g->coeffs))
        fmpz_mod_mpoly_set_coeff_ui_ui(U, 1, expg, ctx);
    else
    {
        fmpz_t d;
        fmpz_init(d);
        fmpz_mod_inv(d, g->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_set_coeff_fmpz_ui(U, d, expg, ctx);
        fmpz_clear(d);
    }

    fmpz_mod_mpoly_mul(U, U, g, ctx);

    fmpz_mod_mpoly_sub(res, T, U, ctx);

    flint_free(exp);
    flint_free(expf);
    flint_free(expg);
    fmpz_mod_mpoly_clear(T, ctx);
    fmpz_mod_mpoly_clear(U, ctx);
}
