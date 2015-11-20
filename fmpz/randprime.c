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

void fmpz_randprime(fmpz_t f, flint_rand_t state, mp_bitcnt_t bits, int proved)
{
    if (bits <= FLINT_BITS-3)
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
            mpz_urandomb(mpz_ptr, state->gmp_state, bits-1);
            mpz_setbit(mpz_ptr, bits-1);
            _fmpz_demote_val(f);
            
            fmpz_nextprime(f, f, proved);
        } while (fmpz_bits(f) != bits);
    }
}
