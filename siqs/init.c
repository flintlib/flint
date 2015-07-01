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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

void qsieve_init(qs_t qs_inf, fmpz_t n)
{
    ulong i;

    /* store n in struct */
    fmpz_init_set(qs_inf->n, n);

    /* determine the number of bits of n */
    qs_inf->bits =  fmpz_bits(n);

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < QS_TUNE_SIZE; i++)
    {
        if (qsieve_tune[i][0] > qs_inf->bits)
            break;
    }
    i--;

    qs_inf->ks_primes  = qsieve_tune[i][1]; /* number of Knuth-Schroeppel primes */
    qs_inf->num_primes = qsieve_tune[i][2]; /* number of factor base primes */

    fmpz_init(qs_inf->kn); /* initialise kn */

    qs_inf->factor_base = NULL;
    qs_inf->sqrts       = NULL;

}
