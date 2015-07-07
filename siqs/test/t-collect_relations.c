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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

int main(void)
{
    slong i, j, ncols, nrows;
    mp_limb_t small_factor;
    char * sieve;
    fmpz_t n, temp, temp2;
    qs_t qs_inf;
    fmpz_init(n);

    FLINT_TEST_INIT(state);

    flint_printf("collect_relations....");
    fflush(stdout);

    for (i = 0; i < 1000; i++)
    {
        fmpz_randtest_unsigned(n, state, 130);

        if (fmpz_is_zero(n) || fmpz_is_one(n) || fmpz_bits(n) <= 60 || fmpz_bits(n) >= 100) continue;

        qsieve_init(qs_inf, n);

        small_factor = qsieve_knuth_schroeppel(qs_inf);

        if (small_factor) continue;

        fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k); /* haven't calculated earlier */

        small_factor = qsieve_primes_init(qs_inf);

        if (small_factor) continue;

        sieve = (char *) flint_malloc(qs_inf->sieve_size + sizeof(ulong));

        qsieve_linalg_init(qs_inf);    /* initialize linear algebra fields */

        qs_inf->sieve_bits = 30; /* sieve threshold */

        qsieve_collect_relations(qs_inf, sieve);  /* perform sieving */

        flint_free(sieve);

        ncols = 0;

        for (j = 0; j < qs_inf->buffer_size; j++)
        {
            if (qs_inf->matrix[j].weight > 0) ncols++;
        }

        flint_printf("total coulmns = %wd\n", ncols);

        ncols = qs_inf->num_primes + qs_inf->extra_rels;
        nrows = qs_inf->num_primes;

        reduce_matrix(qs_inf, &nrows, &ncols, qs_inf->matrix);

        qsieve_linalg_clear(qs_inf);        /* free linear algebra field */
        qsieve_poly_clear(qs_inf);          /* free polynomials field */
        qsieve_clear(qs_inf);               /* free rest of the field */
        fmpz_clear(n);

        break;
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");

    return 0;
}


