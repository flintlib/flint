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
#include "fmpz_mod_mpoly.h"

int
fmpz_mod_mpoly_vec_is_autoreduced(const fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j, len, alloc;
    fmpz_mod_mpoly_t h;
    fmpz_t scale;
    fmpz_mod_mpoly_struct ** Q, ** B;
    int result;

    len = G->length;
    alloc = len - 1;

    if (len == 0)
        return 1;

    fmpz_init(scale);
    Q = flint_malloc(sizeof(fmpz_mod_mpoly_struct *) * alloc);
    B = flint_malloc(sizeof(fmpz_mod_mpoly_struct *) * alloc);

    for (i = 0; i < alloc; i++)
    {
        Q[i] = flint_malloc(sizeof(fmpz_mod_mpoly_struct));
        fmpz_mod_mpoly_init(Q[i], ctx);
    }

    fmpz_mod_mpoly_init(h, ctx);

    result = 1;
    for (i = 0; i < len && result; i++)
    {
        for (j = 0; j < i; j++)
            B[j] = fmpz_mod_mpoly_vec_entry(G, j);
        for (j = i + 1; j < G->length; j++)
            B[j - 1] = fmpz_mod_mpoly_vec_entry(G, j);

        fmpz_mod_mpoly_divrem_ideal(Q, h, fmpz_mod_mpoly_vec_entry(G, i), B, G->length - 1, ctx);
        fmpz_mod_mpoly_make_monic(h, h, ctx);

        if (fmpz_mod_mpoly_is_zero(h, ctx) || !fmpz_mod_mpoly_equal(h, fmpz_mod_mpoly_vec_entry(G, i), ctx))
            result = 0;
    }

    for (i = 0; i < alloc; i++)
    {
        fmpz_mod_mpoly_clear(Q[i], ctx);
        flint_free(Q[i]);
    }

    fmpz_mod_mpoly_clear(h, ctx);

    flint_free(Q);
    flint_free(B);
    fmpz_clear(scale);

    return result;
}
