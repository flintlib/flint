/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "utils_flint.h"

void
fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_t res, slong var, const fmpz_poly_t pol, const fmpz_mpoly_ctx_t ctx)
{
    if (ctx->minfo->nvars == 0)
    {
        flint_printf("fmpz_mpoly_set_gen_fmpz_poly: require nvars >= 1");
        flint_abort();
    }

    if (fmpz_poly_is_zero(pol))
    {
        fmpz_mpoly_zero(res, ctx);
    }
    else if (pol->length == 1)
    {
        fmpz_mpoly_set_fmpz(res, pol->coeffs, ctx);
    }
    else if (pol->length == 2)
    {
        fmpz_mpoly_gen(res, var, ctx);
        fmpz_mpoly_scalar_mul_fmpz(res, res, pol->coeffs + 1, ctx);
        fmpz_mpoly_add_fmpz(res, res, pol->coeffs + 0, ctx);
    }
    else
    {
        slong len, num, i;
        ulong * exp;

        len = pol->length;
        exp = flint_malloc(fmpz_mpoly_ctx_nvars(ctx) * sizeof(ulong)); /* TMP_ALLOC? */

        for (i = 0; i < fmpz_mpoly_ctx_nvars(ctx); i++)
            exp[i] = 0;

        num = 1;
        for (i = pol->length - 2; i >= 0; i--)
            num += !fmpz_is_zero(pol->coeffs + i);

        fmpz_mpoly_fit_bits(res, FLINT_BIT_COUNT(len), ctx);
        fmpz_mpoly_fit_length(res, num, ctx);

        /* xxx? zero fmpzs? */
        res->length = 0;

        /* xxx: reverse ordering */
        for (i = len - 1; i >= 0; i--)
        {
            if (!fmpz_is_zero(pol->coeffs + i))
            {
                exp[var] = i;
                fmpz_mpoly_push_term_fmpz_ui(res, pol->coeffs + i, exp, ctx);
            }
        }

        _fmpz_mpoly_set_length(res, num, ctx);

        flint_free(exp);
    }
}
