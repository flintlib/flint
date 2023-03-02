/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz

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
fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m)
{
    flint_bitcnt_t bits = fmpz_bits(m);
    int sgn = fmpz_sgn(m);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        _fmpz_demote(f);
        *f =  (sgn >= 0) ? n_randint(state, *m) : - n_randint(state, -(*m));
    }
    else
    {
        __mpz_struct * mf = _fmpz_promote(f);

        _flint_rand_init_gmp(state);
        mpz_urandomm(mf, state->gmp_state, COEFF_TO_PTR(*m));
        if (sgn < 0)
            mpz_neg(mf, mf);
        _fmpz_demote_val(f);
    }
}
