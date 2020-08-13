/*
    Authored 2015 by Daniel S. Roche; US Government work in the public domain.

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void fmpz_randprime(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits, int proved)
{
    if (bits <= FLINT_BITS - 2)
    {
        _fmpz_demote(f);
        *f = n_randprime(state, bits, proved);
    }
    else
    {
        /* Here I would like to just call
         * fmpz_randbits(f, state, bits);
         * but it has different semantics from n_randbits,
         * and in particular may return integers with fewer bits.
         */
        __mpz_struct *mpz_ptr = _fmpz_promote(f);
        _flint_rand_init_gmp(state);

        do
        {
            mpz_urandomb(mpz_ptr, state->gmp_state, bits - 1);
            mpz_setbit(mpz_ptr, bits - 1);

            fmpz_nextprime(f, f, proved);
        } while (fmpz_bits(f) != bits);
    }
}
