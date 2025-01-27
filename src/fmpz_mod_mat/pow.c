/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

void
fmpz_mod_mat_pow(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong exp,
                 const fmpz_mod_ctx_t ctx)
{
    slong d = fmpz_mod_mat_nrows(A, ctx);

    if (exp <= 2 || d <= 1)
    {
        if (exp == 0 || d == 0)
        {
            fmpz_mod_mat_one(B, ctx);
        }
        else if (d == 1)
        {
            fmpz_mod_pow_ui(fmpz_mod_mat_entry(B, 0, 0),
                            fmpz_mod_mat_entry(A, 0, 0), exp, ctx);
        }
        else if (exp == 1)
        {
            fmpz_mod_mat_set(B, A, ctx);
        }
        else if (exp == 2)
        {
            fmpz_mod_mat_sqr(B, A, ctx);
        }
    }
    else
    {
        fmpz_mod_mat_t T, U;
        slong i;

        fmpz_mod_mat_init_set(T, A, ctx);
        fmpz_mod_mat_init(U, d, d, ctx);

        for (i = ((slong) FLINT_BIT_COUNT(exp)) - 2; i >= 0; i--)
        {
            fmpz_mod_mat_sqr(U, T, ctx);

            if (exp & (WORD(1) << i))
                fmpz_mod_mat_mul(T, U, A, ctx);
            else
                fmpz_mod_mat_swap(T, U, ctx);
        }

        fmpz_mod_mat_swap(B, T, ctx);
        fmpz_mod_mat_clear(T, ctx);
        fmpz_mod_mat_clear(U, ctx);
    }
}
