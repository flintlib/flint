/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_mul(
    fmpz * A,
    const fmpz * B,
    const fmpz * C,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    for (len--; len >= 0; len--)
        fmpz_mod_mul(A + len, B + len, C + len, ctx);
}

