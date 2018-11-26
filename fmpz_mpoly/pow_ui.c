/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_pow_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong k, const fmpz_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fmpz_mpoly_set_ui(A, k == 0, ctx);
        return;
    }
    else if (k <= 2)
    {
        if (k == 2)
        {
            fmpz_mpoly_mul(A, B, B, ctx);
        }
        else if (k == 1)
        {
            fmpz_mpoly_set(A, B, ctx);
        }
        else
        {
            fmpz_mpoly_set_ui(A, 1, ctx);
        }
    }
    else
    {
        fmpz_mpoly_pow_fps(A, B, k, ctx);
    }
}
