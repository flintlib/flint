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
fmpz_mpoly_vec_set_primitive_unique(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, len;

    fmpz_mpoly_vec_set(G, F, ctx);

    len = G->length;

    for (i = 0; i < len; i++)
    {
        /* skip zero */
        if (fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(G, i), ctx))
        {
            fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, i),
                fmpz_mpoly_vec_entry(G, len - 1), ctx);
            G->length--;
            len--;
            i--;
        }
        else
        {
            fmpz_mpoly_primitive_part(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, i), ctx);

            for (j = 0; j < i; j++)
            {
                if (fmpz_mpoly_equal(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, j), ctx))
                {
                    fmpz_mpoly_zero(fmpz_mpoly_vec_entry(G, i), ctx);
                    fmpz_mpoly_swap(fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, len - 1), ctx);
                    G->length--;
                    len--;
                    i--;
                    break;
                }
            }
        }
    }
}
