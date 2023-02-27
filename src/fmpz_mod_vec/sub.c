/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_sub(fmpz * a, const fmpz * b, const fmpz * c,
                                             slong n, const fmpz_mod_ctx_t ctx)
{
    while (--n >= 0)
        fmpz_mod_sub(a + n, b + n, c + n, ctx);
}
