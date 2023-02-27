/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


static void _fmpz_mpoly_factor_mul_mpoly_fmpz(
    fmpz_mpoly_factor_t f,
    fmpz_mpoly_t A,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mpoly_get_fmpz(t, A, ctx);
        fmpz_pow_fmpz(t, t, e);
        fmpz_mul(f->constant, f->constant, t);
        fmpz_clear(t);
    }
    else
    {
        fmpz_mpoly_factor_append_fmpz_swap(f, A, e, ctx);
    }
}

/*
    append factor_squarefree(A)^e to f
    assuming is A is primitive wrt to each variable
*/
int _fmpz_mpoly_factor_squarefree(
    fmpz_mpoly_factor_t f,
    fmpz_mpoly_t A,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong v;
    fmpz_t k, ke;
    fmpz_mpoly_t S, Sp, Sm, Ss, Y, Z;

    if (A->length < 2)
    {
        _fmpz_mpoly_factor_mul_mpoly_fmpz(f, A, e, ctx);
        return 1;
    }

    fmpz_init(k);
    fmpz_init(ke);
    fmpz_mpoly_init(S, ctx);
    fmpz_mpoly_init(Sp, ctx);
    fmpz_mpoly_init(Sm, ctx);
    fmpz_mpoly_init(Ss, ctx);
    fmpz_mpoly_init(Y, ctx);
    fmpz_mpoly_init(Z, ctx);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_mpoly_derivative(Sp, A, v, ctx);

        if (fmpz_mpoly_is_zero(Sp, ctx))
            continue;

        success = fmpz_mpoly_gcd_cofactors(Sm, Ss, Y, A, Sp, ctx);
        if (!success)
            continue;

        for (fmpz_set_ui(k, 1); !(fmpz_mpoly_derivative(Sp, Ss, v, ctx),
                                  fmpz_mpoly_sub(Z, Y, Sp, ctx),
                                  fmpz_mpoly_is_zero(Z, ctx));
                                                     fmpz_add_ui(k, k, 1))
        {
            success = fmpz_mpoly_gcd_cofactors(S, Ss, Y, Ss, Z, ctx);
            if (!success)
                continue;

            fmpz_mul(ke, k, e);
            _fmpz_mpoly_factor_mul_mpoly_fmpz(f, S, k, ctx);
        }

        fmpz_mul(ke, k, e);
        _fmpz_mpoly_factor_mul_mpoly_fmpz(f, Ss, k, ctx);

        success = 1;
        goto cleanup;
    }

    success = 0;

cleanup:

    fmpz_clear(k);
    fmpz_mpoly_clear(S, ctx);
    fmpz_mpoly_clear(Sp, ctx);
    fmpz_mpoly_clear(Sm, ctx);
    fmpz_mpoly_clear(Ss, ctx);
    fmpz_mpoly_clear(Y, ctx);
    fmpz_mpoly_clear(Z, ctx);

    return success;
}


int fmpz_mpoly_factor_squarefree(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz_mpoly_factor_t g;

    fmpz_mpoly_factor_init(g, ctx);

    success = fmpz_mpoly_factor_content(g, A, ctx);
    if (!success)
        goto cleanup;

    fmpz_swap(f->constant, g->constant);
    f->num = 0;
    for (i = 0; i < g->num; i++)
    {
        success = _fmpz_mpoly_factor_squarefree(f, g->poly + i, g->exp + i, ctx);
        if (!success)
            goto cleanup;
    }

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}

