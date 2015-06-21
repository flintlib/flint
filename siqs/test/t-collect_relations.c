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
    slong i, j;
    mp_limb_t small_factor;
    fmpz_t n, temp, temp2;
    qs_t qs_inf;
    fmpz_init(n);
    fmpz_init(temp);
    fmpz_init(temp2);

    FLINT_TEST_INIT(state);

    flint_printf("collect_relations....");
    fflush(stdout);

    for (i = 0; i < 1000; i++)
    {
        fmpz_randtest_unsigned(n, state, 130);

        if (fmpz_is_zero(n) || fmpz_is_one(n) || fmpz_bits(n) <= 100) continue;

        qsieve_init(qs_inf, n);
        small_factor = qsieve_knuth_schroeppel(qs_inf);

        if (small_factor) continue;

        fmpz_mul_ui(qs_inf->kn, qs_inf->n, qs_inf->k); /* haven't calculated earlier */
        small_factor = qsieve_primes_init(qs_inf);

        if (small_factor) continue;

        flint_printf("number to factor : ");
        fmpz_print(qs_inf->kn);
        flint_printf("\n");

        flint_printf("optimal hypercube: ");
        fmpz_print(qs_inf->target_A);
        flint_printf("\n");

        qsieve_poly_init(qs_inf);

        qsieve_compute_q0(qs_inf);

        qsieve_init_A0(qs_inf);

        qsieve_compute_pre_data(qs_inf);

        qsieve_init_poly_first(qs_inf);

        flint_printf("approximated hypercube: ");
        fmpz_print(qs_inf->A);
        flint_printf("\n with factor = %wu \n", qs_inf->s);

        qsieve_compute_C(qs_inf);

        char * sieve = flint_malloc(qs_inf->sieve_size + sizeof(ulong));

        qsieve_linalg_init(qs_inf);

        qsieve_do_sieving(qs_inf, sieve);

        qs_inf->sieve_bits = 32;

        qsieve_evaluate_sieve(qs_inf, sieve);

        qsieve_next_A0(qs_inf);

        qsieve_compute_pre_data(qs_inf);

        qsieve_init_poly_first(qs_inf);

        flint_printf("approximated hypercube: ");
        fmpz_print(qs_inf->A);
        flint_printf("\n with factor = %wu \n", qs_inf->s);

        qsieve_clear(qs_inf);

        break;
    }

    fmpz_clear(n);
    fmpz_clear(temp);
    fmpz_clear(temp2);

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");

    return 0;
}


