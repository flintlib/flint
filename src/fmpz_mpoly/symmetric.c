/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void
fmpz_mpoly_symmetric_gens(fmpz_mpoly_t res, ulong k, slong * vars, slong n, const fmpz_mpoly_ctx_t ctx)
{
    ulong * exp;
    slong * c;
    slong nvars, i, j;

    if (k == 0)
    {
        fmpz_mpoly_one(res, ctx);
        return;
    }

    fmpz_mpoly_zero(res, ctx);

    if (k > n)
        return;

    nvars = ctx->minfo->nvars;

    /* Generate all combinations (Knuth algorithm L). Todo:
    Knuth algorithm T is faster. */

    c = flint_malloc(sizeof(slong) * (k + 2));
    exp = flint_calloc(nvars, sizeof(ulong));

    for (j = 0; j < k; j++)
        c[j] = j;

    c[k] = n;
    c[k + 1] = 0;

    while (1)
    {
        /* Visit this combination */
        for (i = 0; i < n; i++)
            exp[vars[i]] = 0;

        for (i = 0; i < k; i++)
            exp[vars[c[i]]] = 1;

        fmpz_mpoly_push_term_ui_ui(res, 1, exp, ctx);

        j = 1;
        while (c[j - 1] + 1 == c[j])
        {
            c[j - 1] = j - 1;
            j++;
            if (c[j - 1] + 1 != c[j])
                break;
        }

        if (j > k)
            break;

        c[j - 1]++;
    }

    fmpz_mpoly_sort_terms(res, ctx);

    flint_free(c);
    flint_free(exp);
}

void
fmpz_mpoly_symmetric(fmpz_mpoly_t res, ulong k, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars;
    slong * vars;

    nvars = ctx->minfo->nvars;

    vars = flint_malloc(sizeof(slong) * nvars);
    for (i = 0; i < nvars; i++)
        vars[i] = i;
    fmpz_mpoly_symmetric_gens(res, k, vars, nvars, ctx);
    flint_free(vars);
}
