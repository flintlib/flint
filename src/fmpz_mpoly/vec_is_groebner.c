/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int
fmpz_mpoly_vec_is_groebner(const fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, len;
    fmpz_mpoly_t h;
    int result;

    len = G->length;

    if (len == 0)
        return 1;

    fmpz_mpoly_init(h, ctx);
    result = 1;

    for (i = 0; i < len && result; i++)
    {
        for (j = i + 1; j < len && result; j++)
        {
            fmpz_mpoly_spoly(h, fmpz_mpoly_vec_entry(G, i), fmpz_mpoly_vec_entry(G, j), ctx);

            fmpz_mpoly_reduction_primitive_part(h, h, G, ctx);
            if (!fmpz_mpoly_is_zero(h, ctx))
                result = 0;
        }
    }

    if (F != NULL)
    {
        for (i = 0; i < F->length && result; i++)
        {
            fmpz_mpoly_reduction_primitive_part(h, fmpz_mpoly_vec_entry(F, i), G, ctx);
            if (!fmpz_mpoly_is_zero(h, ctx))
                result = 0;
        }
    }

    fmpz_mpoly_clear(h, ctx);
    return result;
}
