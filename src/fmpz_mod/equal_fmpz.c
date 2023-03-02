/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int fmpz_mod_equal_fmpz(const fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
{
    int res;
    fmpz_t t;
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
    fmpz_init(t);
    fmpz_mod_set_fmpz(t, b, ctx);
    res = fmpz_equal(a, t);
    fmpz_clear(t);
    return res;
}

int fmpz_mod_equal_si(const fmpz_t a, slong b, const fmpz_mod_ctx_t ctx)
{
    int res;
    fmpz_t t;
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
    fmpz_init(t);
    fmpz_mod_set_si(t, b, ctx);
    res = fmpz_equal(a, t);
    fmpz_clear(t);
    return res;
}
