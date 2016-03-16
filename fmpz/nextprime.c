/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

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
    else if (FLINT_BIT_COUNT(*n) < FLINT_BITS - 2)
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
        __mpz_struct *res_mpz = _fmpz_promote(res);
        mpz_nextprime(res_mpz, res_mpz);
        _fmpz_demote_val(res);
    }

    if (proved)
    {
        int proof = fmpz_is_prime(res);
        
        if (proof == 0)
        {
            /* Keep searching. No big penalty for recursion here because this
             * will almost never happen.
             */
            fmpz_nextprime(res, res, proved);
        }
        else if (proof < 0)
        {
            flint_printf("Exception in fmpz_nextprime: Proof requested but couldn't be found\n");
            flint_abort();
        }
    }
}

