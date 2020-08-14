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

int fmpz_mpoly_factor_squarefree(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, v;
    fmpz_mpoly_t c;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_factor_t g;
    fmpz * var_powers;
    fmpz_t k;
    fmpz_mpoly_t S, Sp, Sm, Ss, Y, Z;

    f->num = 0;
    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_mpoly_get_fmpz(f->constant, A, ctx);
        return 1;
    }

    FLINT_ASSERT(A->length > 0);
    _fmpz_vec_content(f->constant, A->coeffs, A->length);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_neg(f->constant, f->constant);

    fmpz_mpoly_factor_fit_length(f, 1, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(f->poly + 0, A, f->constant, ctx);
    fmpz_one(f->exp + 0);
    f->num = 1;

    fmpz_init(k);
    fmpz_mpoly_init(S, ctx);
    fmpz_mpoly_init(Sp, ctx);
    fmpz_mpoly_init(Sm, ctx);
    fmpz_mpoly_init(Ss, ctx);
    fmpz_mpoly_init(Y, ctx);
    fmpz_mpoly_init(Z, ctx);

    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(c, ctx);
    var_powers = _fmpz_vec_init(ctx->minfo->nvars);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        fmpz_set(g->constant, f->constant);
        g->num = 0;
        for (j = 0; j < f->num; j++)
        {
            FLINT_ASSERT(fmpz_is_one(f->exp + j));

            fmpz_mpoly_to_univar(u, f->poly + j, v, ctx);
            FLINT_ASSERT(u->length > 0);

            success = _fmpz_mpoly_vec_content_mpoly(c, u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            fmpz_add(var_powers + v, var_powers + v, u->exps + u->length - 1);
            for (i = 0; i < u->length; i++)
            {
                fmpz_sub(u->exps + i, u->exps + i, u->exps + u->length - 1);
                success = fmpz_mpoly_divides(u->coeffs + i, u->coeffs + i, c, ctx);
                FLINT_ASSERT(success);
            }

            fmpz_one(k);
            _fmpz_mpoly_factor_mul_mpoly_fmpz(g, c, k, ctx);

            if (u->length > 1)
            {
                fmpz_mpoly_from_univar_bits(c, A->bits, u, v, ctx);
                fmpz_mpoly_factor_append_ui(g, c, 1, ctx);
            }
            else
            {
                FLINT_ASSERT(fmpz_mpoly_is_one(u->coeffs + 0, ctx));
            }
        }

        fmpz_mpoly_factor_swap(f, g, ctx);
    }

    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        int loop_failed = 0;

        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            fmpz_mpoly_derivative(Sp, f->poly + j, v, ctx);

            if (fmpz_mpoly_is_zero(Sp, ctx))
                continue;

            success = fmpz_mpoly_gcd_cofactors(Sm, Ss, Y, f->poly + j, Sp, ctx);
            if (!success)
            {
                loop_failed = 1;
                continue;
            }

            for (fmpz_set_ui(k, 1); !(fmpz_mpoly_derivative(Sp, Ss, v, ctx),
                                      fmpz_mpoly_sub(Z, Y, Sp, ctx),
                                      fmpz_mpoly_is_zero(Z, ctx));
                                                         fmpz_add_ui(k, k, 1))
            {
                success = fmpz_mpoly_gcd_cofactors(S, Ss, Y, Ss, Z, ctx);
                if (!success)
                    goto cleanup;

                _fmpz_mpoly_factor_mul_mpoly_fmpz(g, S, k, ctx);
            }

            _fmpz_mpoly_factor_mul_mpoly_fmpz(g, Ss, k, ctx);
            goto continue_outer;
        }

        if (loop_failed)
        {
            success = 0;
            goto cleanup;
        }

        FLINT_ASSERT(fmpz_mpoly_is_fmpz(f->poly + j, ctx));
        fmpz_one(k);
        _fmpz_mpoly_factor_mul_mpoly_fmpz(f, f->poly + j, k, ctx);

continue_outer:
        (void) NULL;
    }

    fmpz_mpoly_factor_swap(f, g, ctx);

    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        if (fmpz_is_zero(var_powers + v))
            continue;

        fmpz_mpoly_factor_fit_length(f, f->num + 1, ctx);
        fmpz_mpoly_gen(f->poly + f->num, v, ctx);
        fmpz_swap(f->exp + f->num, var_powers + v);
        f->num++;
    }

    success = 1;

cleanup:

    fmpz_clear(k);
    fmpz_mpoly_clear(S, ctx);
    fmpz_mpoly_clear(Sp, ctx);
    fmpz_mpoly_clear(Sm, ctx);
    fmpz_mpoly_clear(Ss, ctx);
    fmpz_mpoly_clear(Y, ctx);
    fmpz_mpoly_clear(Z, ctx);

    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(c, ctx);
    _fmpz_vec_clear(var_powers, ctx->minfo->nvars);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}
