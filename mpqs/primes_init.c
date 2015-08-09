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

    Built upon existing FLINT siqs
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

mp_limb_t
mpqs_primes_init(mpqs_t mpqs_inf)
{
    slong num_primes;
    slong i;
    mp_limb_t k = mpqs_inf->k;
    mp_limb_t small_factor = 0;
    prime_t * factor_base;

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < MPQS_TUNE_SIZE; i++)
    {
        if (mpqs_tune[i][0] > mpqs_inf->bits)
            break;
    }
    i--;

    num_primes = mpqs_tune[i][2]; /* number of factor base primes */
    mpqs_inf->sieve_size = mpqs_tune[i][4]; /* size of sieve to use */
    mpqs_inf->small_primes = mpqs_tune[i][3]; /* number of primes to not sieve with */

    mpqs_inf->num_primes = 0; /* start with 0 primes */
    factor_base = mpqs_compute_factor_base(&small_factor, mpqs_inf, num_primes + mpqs_inf->ks_primes); /* build up FB */
    
    if (small_factor)
        return small_factor;

    mpqs_inf->num_primes = num_primes;

    /* consider k and 2 as factor base primes */
    factor_base[0].p = k;
    factor_base[0].pinv = n_preinvert_limb(k);
    factor_base[0].size = FLINT_BIT_COUNT(k);
    factor_base[1].p = 2;
    factor_base[0].size = 2;

    return 0;
}
