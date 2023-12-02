/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <limits.h>

#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

void
fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    ulong m;

    fmpz_randtest_unsigned(f, state, bits);

    m = n_randlimb(state);
    if (m & UWORD(1))
        fmpz_neg(f, f);
}

void
fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    ulong m;

    m    = n_randlimb(state);
    bits = n_randint(state, bits + 1);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        _fmpz_demote(f);
        if (m & UWORD(3))
            *f = n_randtest_bits(state, bits);
        else
        {
            m >>= 2;
            if (bits == 0)
                *f = 0;
            else if (bits < SMALL_FMPZ_BITCOUNT_MAX)
                *f = m & UWORD(1);
            else
                *f = COEFF_MAX;
        }
    }
    else
    {
        __mpz_struct * mf = _fmpz_promote(f);

        _flint_rand_init_gmp(state);
        mpz_rrandomb(mf, state->gmp_state, bits);
        _fmpz_demote_val(f);
    }
}

void
fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    if (bits == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_randtest_not_zero). bits == 0.\n");
    }

    fmpz_randtest(f, state, bits);
    if (fmpz_is_zero(f))
        fmpz_one(f);
}
