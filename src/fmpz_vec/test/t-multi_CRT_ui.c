/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "nmod_vec.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_multi_CRT_ui, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong i, j, len, pbits, num_primes;
        ulong p;
        fmpz_t M;
        fmpz * A, * B;
        int sign;
        nn_ptr primes, Amod;
        nn_srcptr * Amod_ptr;

        flint_set_num_threads(1 + n_randint(state, 4));

        if (n_randint(state, 2) == 0)
        {
            len = n_randint(state, 200);
            num_primes = n_randint(state, 40);
            pbits = 1;
        }
        else
        {
            len = n_randint(state, 10);
            num_primes = n_randint(state, 4);
            pbits = n_randint(state, FLINT_BITS);
        }

        sign = n_randint(state, 2);

        A = _fmpz_vec_init(len);
        B = _fmpz_vec_init(len);
        primes = _nmod_vec_init(num_primes);
        Amod = _nmod_vec_init(num_primes * len);
        Amod_ptr = flint_malloc(sizeof(nn_ptr) * num_primes);

        for (i = 0; i < num_primes; i++)
            Amod_ptr[i] = Amod + i * len;

        fmpz_init(M);
        fmpz_one(M);

        p = n_nextprime(UWORD(1) << pbits, 0);

        for (i = 0; i < num_primes; i++)
        {
            primes[i] = p;
            fmpz_mul_ui(M, M, p);
            p = n_nextprime(p, 0);
        }

        for (i = 0; i < len; i++)
        {
            if (sign && !fmpz_is_one(M))
                fmpz_randtest_mod_signed(A + i, state, M);
            else
                fmpz_randtest_mod(A + i, state, M);

            for (j = 0; j < num_primes; j++)
                Amod[j * len + i] = fmpz_fdiv_ui(A + i, primes[j]);
        }

        _fmpz_vec_multi_CRT_ui(B, Amod_ptr, len, primes, num_primes, sign);

        if (!_fmpz_vec_equal(A, B, len))
        {
            flint_printf("FAIL\n");
            flint_printf("primes = %{ulong*}\n\n", primes, num_primes);
            flint_printf("A = %{fmpz*}\n\n", A, len);
            flint_printf("B = %{fmpz*}\n\n", B, len);
            flint_abort();
        }

        _fmpz_vec_clear(A, len);
        _fmpz_vec_clear(B, len);
        _nmod_vec_clear(primes);
        _nmod_vec_clear(Amod);
        flint_free(Amod_ptr);
        fmpz_clear(M);
    }

    TEST_FUNCTION_END(state);
}
