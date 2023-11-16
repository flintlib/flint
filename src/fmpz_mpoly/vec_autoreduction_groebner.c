/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

static int
monomial_divides(const ulong * exp1, const ulong * exp2, slong nvars)
{
    slong i;
    for (i = 0; i < nvars; i++)
        if (exp1[i] > exp2[i])
            return 0;
    return 1;
}

void
fmpz_mpoly_vec_autoreduction_groebner(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t Gin, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, nvars;
    ulong * exp1, * exp2;

    if (G != Gin)
        fmpz_mpoly_vec_set(G, Gin, ctx);

    /* Content removal */
    for (i = 0; i < G->length; i++)
        fmpz_mpoly_primitive_part(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, i), ctx);

    /* Make sure there are no duplicate or zero entries */
    for (i = 0; i < G->length; i++)
    {
        if (fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(G, i), ctx))
        {
            fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, G->length - 1), ctx);
            fmpz_mpoly_vec_set_length(G, G->length - 1, ctx);
        }
        else
        {
            for (j = i + 1; j < G->length; j++)
            {
                if (fmpz_mpoly_equal(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, j), ctx))
                {
                    fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, j), fmpz_mpoly_vec_entry(G, G->length - 1), ctx);
                    fmpz_mpoly_vec_set_length(G, G->length - 1, ctx);
                }
            }
        }
    }

    if (G->length <= 1)
        return;

    /* First filter based on divisibility of leading terms */

    nvars = ctx->minfo->nvars;
    exp1 = flint_malloc(nvars * sizeof(ulong));
    exp2 = flint_malloc(nvars * sizeof(ulong));

    for (i = 0; i < G->length; i++)
    {
        fmpz_mpoly_get_term_exp_ui(exp1, fmpz_mpoly_vec_entry(G, i), 0, ctx);

        for (j = 0; j < G->length; j++)
        {
            if (i != j)
            {
                fmpz_mpoly_get_term_exp_ui(exp2, fmpz_mpoly_vec_entry(G, j), 0, ctx);

                if (monomial_divides(exp2, exp1, nvars))
                {
                    fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, G->length - 1), ctx);
                    fmpz_mpoly_vec_set_length(G, G->length - 1, ctx);
                    break;
                }
            }
        }
    }

    flint_free(exp1);
    flint_free(exp2);

    /* Now do inter-reduction */
    if (G->length >= 2)
    {
        fmpz_t scale;
        fmpz_mpoly_struct ** Q, ** B;
        slong i, j, alloc;

        alloc = G->length - 1;

        fmpz_init(scale);
        Q = flint_malloc(sizeof(fmpz_mpoly_struct *) * alloc);
        B = flint_malloc(sizeof(fmpz_mpoly_struct *) * alloc);

        for (i = 0; i < alloc; i++)
        {
            Q[i] = flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(Q[i], ctx);
        }

        for (i = 0; i < G->length; i++)
        {
            for (j = 0; j < i; j++)
                B[j] = fmpz_mpoly_vec_entry(G, j);
            for (j = i + 1; j < G->length; j++)
                B[j - 1] = fmpz_mpoly_vec_entry(G, j);

            fmpz_mpoly_quasidivrem_ideal(scale, Q, fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, i), B, G->length - 1, ctx);
            fmpz_mpoly_primitive_part(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, i), ctx);

            if (fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(G, i), ctx))
            {
                fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, G->length - 1), ctx);
                fmpz_mpoly_vec_set_length(G, G->length - 1, ctx);
                i--;
            }
        }

        for (i = 0; i < alloc; i++)
        {
            fmpz_mpoly_clear(Q[i], ctx);
            flint_free(Q[i]);
        }

        flint_free(Q);
        flint_free(B);
        fmpz_clear(scale);
    }
}
