/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_ctx_init_rand_bits(fmpz_mod_ctx_t ctx,
                                   flint_rand_t state, flint_bitcnt_t max_bits)
{
    fmpz_t m;
    fmpz_init(m);
    fmpz_randtest_unsigned(m, state, max_bits);
    fmpz_add_ui(m, m, 2);
    fmpz_mod_ctx_init(ctx, m);
    fmpz_clear(m);
}

void fmpz_mod_ctx_init_rand_bits_prime(fmpz_mod_ctx_t ctx,
                                   flint_rand_t state, flint_bitcnt_t max_bits)
{
    fmpz_t m;
    fmpz_init(m);
    fmpz_randtest_unsigned(m, state, max_bits);
    fmpz_nextprime(m, m, 0);
    fmpz_mod_ctx_init(ctx, m);
    fmpz_clear(m);
}

