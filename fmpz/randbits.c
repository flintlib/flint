/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    if (bits <= FLINT_BITS - 2)
    {
        _fmpz_demote(f);
        *f = n_randbits(state, bits);
        if (n_randint(state, 2))
            *f = -*f;
    }
    else
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        _flint_rand_init_gmp(state);
        mpz_urandomb(mpz_ptr, state->gmp_state, bits);
        mpz_setbit(mpz_ptr, bits - 1);

        if (n_randint(state, 2))
            mpz_neg(mpz_ptr, mpz_ptr);

        _fmpz_demote_val(f);
    }
}
