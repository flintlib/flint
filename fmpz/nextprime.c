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
#include "ulong_extras.h"

void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved)
{
    if (fmpz_sgn(n) <= 0)
    {
        fmpz_set_ui(res, UWORD(2));
        return;
    }
    else if (COEFF_IS_MPZ(*n))
    {
        /* n is big */
        __mpz_struct *res_mpz = _fmpz_promote(res);
        mpz_nextprime(res_mpz, COEFF_TO_PTR(*n));
    }
    else if (FLINT_BIT_COUNT(*n) < SMALL_FMPZ_BITCOUNT_MAX)
    {
        /* n and res will both be small */
        _fmpz_demote(res);
        *res = n_nextprime(*n, proved);
        return;
    }
    else if (res != n)
    {
        /* n is small, but res might not be */
        mpz_t temp_n;
        __mpz_struct *res_mpz = _fmpz_promote(res);
        flint_mpz_init_set_ui(temp_n, *n);
        mpz_nextprime(res_mpz, temp_n);
        _fmpz_demote_val(res);
        mpz_clear(temp_n);
    }
    else
    {
        /* same as above case, but need to handle aliasing here. */
        __mpz_struct *res_mpz = _fmpz_promote_val(res);
        mpz_nextprime(res_mpz, res_mpz);
        _fmpz_demote_val(res);
    }

    if (proved)
    {
        if (!fmpz_is_prime(res))
        {
            /* Keep searching. No big penalty for recursion here because this
             * will almost never happen.
             */
            fmpz_nextprime(res, res, proved);
        }
    }
}
