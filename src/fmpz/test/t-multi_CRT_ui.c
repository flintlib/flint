/*
    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_multi_CRT_ui, state)
{
    fmpz_t input, temp, prod;
    mp_limb_t * output;
    slong i, j, k;
    flint_bitcnt_t bits1, bits2, bits;
    mp_limb_t * primes;
    fmpz * primes2;
    slong num_primes;
    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int sign = n_randint(state, 2);

        num_primes = 1 + n_randint(state, 400);

        fmpz_init(temp);
        fmpz_init(input);
        fmpz_init(prod);
        output = FLINT_ARRAY_ALLOC(num_primes, mp_limb_t);
        primes = FLINT_ARRAY_ALLOC(num_primes, mp_limb_t);
        primes2 = _fmpz_vec_init(num_primes);

try_again:

        bits1 = n_randint(state, FLINT_BITS*3/4) + FLINT_BITS/4;
        bits2 = n_randint(state, FLINT_BITS*3/4) + FLINT_BITS/4;
        if (bits1 > bits2)
            FLINT_SWAP(ulong, bits1, bits2);

        fmpz_one(prod);
        for (j = 0; j < num_primes; j++)
        {
            bits = bits1 + n_randint(state, bits2 - bits1 + 1);
            primes[j] = n_randbits(state, bits);
            primes[j] = n_nextprime(primes[j], 1);
            fmpz_mul_ui(prod, prod, primes[j]);
            fmpz_set_ui(primes2 + j, primes[j]);
        }

        _fmpz_vec_sort(primes2, num_primes);

        for (j = 1; j < num_primes; j++)
            if (fmpz_equal(primes2 + j, primes2 + j - 1))
                goto try_again;

        fmpz_randtest(input, state, n_randint(state, 300) + 1);

        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_temp_init(comb_temp, comb);
        fmpz_multi_mod_ui(output, input, comb, comb_temp);
        fmpz_multi_CRT_ui(temp, output, comb, comb_temp, sign);
        fmpz_comb_temp_clear(comb_temp);
        fmpz_comb_clear(comb);

        if (sign ? fmpz_cmp2abs(prod, temp) < 0 :
                   (fmpz_sgn(temp) < 0 || fmpz_cmp(prod, temp) <= 0))
        {
            flint_printf("FAIL: check crt output range\n");
            flint_printf("i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_sub(temp, temp, input);
        if (!fmpz_divisible(temp, prod))
        {
            flint_printf("FAIL: check crt modulo product of primes\n");
            flint_printf("i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        for (k = 0; k < num_primes; k++)
        {
            if (output[k] != fmpz_fdiv_ui(input, primes[k]))
            {
                flint_printf("FAIL: check multi_mod_ui output");
                flint_printf("i = %wd, k = %wd\n", i, k);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(temp);
        fmpz_clear(input);
        fmpz_clear(prod);
        flint_free(output);
        flint_free(primes);
        _fmpz_vec_clear(primes2, num_primes);
    }

    TEST_FUNCTION_END(state);
}
