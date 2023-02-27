/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_scalar_addmul_fmpz_mod(
    fmpz * A,
    const fmpz * B,
    slong len,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    if (fmpz_is_zero(c))
    {
    }
    else if (fmpz_is_one(c))
    {
        _fmpz_mod_vec_add(A, A, B, len, ctx);
    }
    else
    {
        for (len--; len >= 0; len--)
        {
            fmpz_addmul(A + len, B + len, c);
            fmpz_mod_set_fmpz(A + len, A + len, ctx);
        }
    }
}

