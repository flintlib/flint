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

    Copyright (C) 2007 William Hart and David Harvey
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("CRT_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        long bits, prime_bits, rows, cols, num_primes, j;
        fmpz_t mod;
        fmpz_mat_t A, B, C;
        nmod_mat_t Amod;
        mp_limb_t primes[1000];

        bits = n_randint(state, 500) + 1;
        rows = n_randint(state, 10);
        cols = n_randint(state, 10);
        prime_bits = 1 + n_randint(state, FLINT_BITS - 1);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);

        fmpz_mat_randtest(A, state, bits);

        fmpz_init(mod);
        num_primes = 0;
        primes[0] = n_nextprime(1UL << prime_bits, 0);
        fmpz_set_ui(mod, primes[0]);

        /* + 1 for sign */
        while (fmpz_bits(mod) <= bits + 1)
        {
            primes[num_primes + 1] = n_nextprime(primes[num_primes], 0);
            fmpz_mul_ui(mod, mod, primes[num_primes + 1]);
            num_primes++;
        }

        num_primes++;

        nmod_mat_init(Amod, rows, cols, primes[0]);
        fmpz_mat_get_nmod_mat(Amod, A);
        fmpz_mat_set_nmod_mat(B, Amod);
        fmpz_set_ui(mod, primes[0]);

        for (j = 1; j < num_primes; j++)
        {
            nmod_mat_clear(Amod);
            nmod_mat_init(Amod, rows, cols, primes[j]);
            fmpz_mat_get_nmod_mat(Amod, A);
            fmpz_mat_CRT_ui(B, B, mod, Amod, 1);
            fmpz_mul_ui(mod, mod, primes[j]);
        }

        if (!fmpz_mat_equal(B, A))
        {
            printf("FAIL!\n");
            printf("primes: ");
            for (j = 0; j < num_primes; j++)
                printf("%lu ", primes[j]);
            printf("\nA: \n");
            fmpz_mat_print_pretty(A);
            printf("\nB: \n");
            fmpz_mat_print_pretty(B);
            printf("\n");
            abort();
        }

        nmod_mat_clear(Amod);
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_clear(mod);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
