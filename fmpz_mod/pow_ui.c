/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow, const fmpz_mod_ctx_t ctx)
{
    fmpz_t r, x;

    fmpz_init_set_ui(r, 1);
    fmpz_init_set(x, b);

    while (pow != 0)
    {
        if ((pow & 1) != 0)
        {
            fmpz_mod_mul(r, r, x, ctx);
        }
        pow = pow >> 1;
        if (pow != 0)
        {
            fmpz_mod_mul(x, x, x, ctx);
        }
    }

    fmpz_swap(a, r);
    fmpz_clear(r);
    fmpz_clear(x);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
