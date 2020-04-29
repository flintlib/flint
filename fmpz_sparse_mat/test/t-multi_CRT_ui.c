/*
    Copyright (C) 2007 William Hart and David Harvey
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse_mat.h"
#include "nmod_sparse_mat.h"

int
main(void)
{
    slong rep, r, c, min_nnz, max_nnz, bits, prime_bits, num_primes, j, nreps = 1000;
    mp_limb_t primes[1000];
    nmod_t nmod;
    fmpz_t mod;
    nmod_sparse_mat_struct Amod[1000];
    fmpz_sparse_mat_t A, B;
    FLINT_TEST_INIT(state);

    flint_printf("multi_CRT_ui....");
    fflush(stdout);


    for (rep = 0; rep < nreps; rep++)
    {
        bits = n_randint(state, 500) + 2;
        r = n_randint(state, 10);
        c = n_randint(state, 10);
        min_nnz = 0;
        max_nnz = c;
        prime_bits = 1 + n_randint(state, FLINT_BITS - 1);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(B, r, c);

        fmpz_sparse_mat_randtest(A, state, min_nnz, max_nnz, bits);

        fmpz_init_set_ui(mod, UWORD(1));
        for (num_primes = 0; fmpz_bits(mod) <= bits + 1; num_primes++)
        {
            primes[num_primes] = n_nextprime((num_primes == 0) ? (UWORD(1) << prime_bits) : primes[num_primes - 1], 0);
            nmod_init(&nmod, primes[num_primes]);
            nmod_sparse_mat_init(&Amod[num_primes], r, c, nmod);
            fmpz_sparse_mat_get_nmod_sparse_mat(&Amod[num_primes], A);
            fmpz_mul_ui(mod, mod, primes[num_primes]);
        }
        fmpz_sparse_mat_multi_mod_ui(Amod, num_primes, A);
        fmpz_sparse_mat_multi_CRT_ui(B, Amod, num_primes, 1);

        if (!fmpz_sparse_mat_equal(B, A))
        {
            flint_printf("FAIL!\n");
            flint_printf("primes: ");
            for (j = 0; j < num_primes; j++)
                flint_printf("%wu ", primes[j]);
            flint_printf("\nA: \n");
            fmpz_sparse_mat_print_pretty(A);
            flint_printf("\nB: \n");
            fmpz_sparse_mat_print_pretty(B);
            flint_printf("\n");
            abort();
        }

        for (j = 0; j < num_primes; j++)
            nmod_sparse_mat_clear(&Amod[j]);
        fmpz_sparse_mat_clear(A);
        fmpz_sparse_mat_clear(B);
        fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
