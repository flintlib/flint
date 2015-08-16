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

void
mpqs_init(mpqs_t mpqs_inf, fmpz_t n)
{
    ulong i;

    /* store n in struct */
    fmpz_init_set(mpqs_inf->n, n);

    /* determine the number of bits of n */
    mpqs_inf->bits =  fmpz_bits(n);

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < MPQS_TUNE_SIZE; i++)
    {
        if (mpqs_tune[i][0] > mpqs_inf->bits)
            break;
    }
    i--;

    mpqs_inf->ks_primes  = mpqs_tune[i][1]; /* number of Knuth-Schroeppel primes */
    mpqs_inf->num_primes = mpqs_tune[i][2]; /* number of factor base primes */
    mpqs_inf->qsort_rels = mpqs_tune[i][1];

    fmpz_init(mpqs_inf->kn); /* initialise kn */

    mpqs_inf->factor_base = NULL;
    mpqs_inf->sqrts       = NULL;

}
