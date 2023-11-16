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
fmpz_mpoly_vec_autoreduction(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t Gin, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

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

    /* Now do inter-reduction */
    if (G->length >= 2)
    {
        fmpz_t scale;
        fmpz_mpoly_struct ** Q, ** B;
        fmpz_mpoly_t h;
        slong alloc;
        int changed;

        alloc = G->length - 1;
        fmpz_init(scale);
        fmpz_mpoly_init(h, ctx);
        Q = flint_malloc(sizeof(fmpz_mpoly_struct *) * alloc);
        B = flint_malloc(sizeof(fmpz_mpoly_struct *) * alloc);

        for (i = 0; i < alloc; i++)
        {
            Q[i] = flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(Q[i], ctx);
        }

        while (G->length >= 2)
        {
            changed = 0;

            for (i = 0; i < G->length; i++)
            {
                for (j = 0; j < i; j++)
                    B[j] = fmpz_mpoly_vec_entry(G, j);
                for (j = i + 1; j < G->length; j++)
                    B[j - 1] = fmpz_mpoly_vec_entry(G, j);

                fmpz_mpoly_quasidivrem_ideal(scale, Q, h, fmpz_mpoly_vec_entry(G, i), B, G->length - 1, ctx);
                fmpz_mpoly_primitive_part(h, h, ctx);

                if (!fmpz_mpoly_equal(h, fmpz_mpoly_vec_entry(G, i), ctx))
                {
                    changed = 1;
                    fmpz_mpoly_swap(h, fmpz_mpoly_vec_entry(G, i), ctx);
                }

                if (fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(G, i), ctx))
                {
                    fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, G->length - 1), ctx);
                    fmpz_mpoly_vec_set_length(G, G->length - 1, ctx);
                    i--;
                }
            }

            if (!changed)
                break;
        }

        for (i = 0; i < alloc; i++)
        {
            fmpz_mpoly_clear(Q[i], ctx);
            flint_free(Q[i]);
        }

        flint_free(Q);
        flint_free(B);
        fmpz_clear(scale);
        fmpz_mpoly_clear(h, ctx);
    }
}
