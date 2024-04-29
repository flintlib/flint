/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"

void flint_mpn_rrandom(mp_ptr rp, flint_rand_t state, mp_size_t n)
{
    __mpz_struct str;
    str._mp_d = rp;
    str._mp_alloc = n;
    str._mp_size = n;
    _flint_rand_init_gmp(state);
    FLINT_ASSERT(n >= 1);
    /* Randomly generate numbers where the top limb has the highest bit
       set or not. TODO: do we want this function to generate
       numbers where one or more leading limbs can be zero, or should
       that be a separate function? The documentation does not specify
       precisely what this function does. */
    if (n_randint(state, 2))
        mpz_rrandomb(&str, state->gmp_state, FLINT_BITS * n);
    else
        mpz_rrandomb(&str, state->gmp_state, FLINT_BITS * n - n_randint(state, FLINT_BITS));
}

void flint_mpn_urandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n)
{
    __mpz_struct str;
    str._mp_d = rp;
    str._mp_alloc = (n + FLINT_BITS - 1)/FLINT_BITS;
    str._mp_size = (n + FLINT_BITS - 1)/FLINT_BITS;
    _flint_rand_init_gmp(state);
    mpz_urandomb(&str, state->gmp_state, n);
}