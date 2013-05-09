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

    Copyright (C) 2006, 2011 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

void qsieve_ll_init(qs_t qs_inf, mp_limb_t hi, mp_limb_t lo)
{
    ulong i;
    
    /* store n in struct */
    qs_inf->hi = hi;
    qs_inf->lo = lo;

    /* determine the number of bits of n */
    qs_inf->bits = (hi ? FLINT_BITS + FLINT_BIT_COUNT(hi) : FLINT_BIT_COUNT(lo));

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < QS_LL_TUNE_SIZE; i++)
    {
        if (qsieve_ll_tune[i][0] > qs_inf->bits)
            break;
    }
    i--;
    
    qs_inf->ks_primes  = qsieve_ll_tune[i][1]; /* number of Knuth-Schroeppel primes */
    qs_inf->num_primes = qsieve_ll_tune[i][2]; /* number of factor base primes */

    fmpz_init(qs_inf->kn); /* initialise kn */
    fmpz_init(qs_inf->C); /* initialise C */

    qs_inf->factor_base = NULL;
    qs_inf->sqrts       = NULL;
    qs_inf->B_terms     = NULL;
    qs_inf->A_inv       = NULL;
    qs_inf->A_inv2B     = NULL;

    qs_inf->small       = NULL;
    qs_inf->factor      = NULL;
    qs_inf->matrix      = NULL;
    qs_inf->Y_arr       = NULL;
    qs_inf->relation    = NULL;
    qs_inf->qsort_arr   = NULL;

    qs_inf->prime_count = NULL;

    qs_inf->A = 0;

#if (QS_DEBUG & 16)
    qs_inf->sieve_tally = flint_malloc(256*sizeof(len_t));
#endif
}
