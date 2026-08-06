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

void
fmpz_mod_mpoly_reduction_monic_part(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_vec_t I, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t scale;
    fmpz_mod_mpoly_struct ** Q, ** B;
    slong i, len;

    len = I->length;

    fmpz_init(scale);
    Q = flint_malloc(sizeof(fmpz_mod_mpoly_struct *) * len);
    B = flint_malloc(sizeof(fmpz_mod_mpoly_struct *) * len);

    for (i = 0; i < len; i++)
    {
        Q[i] = flint_malloc(sizeof(fmpz_mod_mpoly_struct));
        fmpz_mod_mpoly_init(Q[i], ctx);
        B[i] = fmpz_mod_mpoly_vec_entry(I, i);
    }

    fmpz_mod_mpoly_divrem_ideal(Q, res, f, B, len, ctx);
    
    if (!fmpz_mod_mpoly_is_zero(res, ctx))
        fmpz_mod_mpoly_make_monic(res, res, ctx);

    for (i = 0; i < len; i++)
    {
        fmpz_mod_mpoly_clear(Q[i], ctx);
        flint_free(Q[i]);
    }

    flint_free(Q);
    flint_free(B);
    fmpz_clear(scale);
}
