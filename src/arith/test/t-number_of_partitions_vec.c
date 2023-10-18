/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"

TEST_FUNCTION_START(arith_number_of_partitions_vec, state)
{
    fmpz * p;
    mp_ptr pmod;
    slong k, n;

    const slong maxn = 1000;


    p = _fmpz_vec_init(maxn);
    pmod = _nmod_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        fmpz_t s, t;
        nmod_t mod;
        nmod_init(&mod, n_randtest_prime(state, 0));

        arith_number_of_partitions_vec(p, n);
        arith_number_of_partitions_nmod_vec(pmod, n, mod);

        for (k = 0; k < n; k++)
        {
            if (fmpz_fdiv_ui(p + k, mod.n) != pmod[k])
            {
                flint_printf("FAIL:\n");
                flint_printf("n = %wd, k = %wd\n", n, k);
                fflush(stdout);
                flint_abort();
            }
        }

        if (n > 1)
        {
            fmpz_init(s);
            fmpz_init(t);

            for (k = 1; k < n; k++)
            {
                slong j;

                j = n - 1 - k*(3*k - 1)/2;
                if (j >= 0)
                    fmpz_set(t, p + j);
                else
                    fmpz_zero(t);

                j = n - 1 - k*(3*k + 1)/2;
                if (j >= 0)
                    fmpz_add(t, t, p + j);

                if (k % 2)
                    fmpz_add(s, s, t);
                else
                    fmpz_sub(s, s, t);
            }

            if (!fmpz_equal(s, p + n - 1))
            {
                flint_printf("FAIL:\n");
                flint_printf("n = %wd\n", n);
                fmpz_print(s);
                flint_printf("\n");
                fmpz_print(p + n - 1);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(s);
            fmpz_clear(t);
        }
    }

    _fmpz_vec_clear(p, maxn);
    _nmod_vec_clear(pmod);

    TEST_FUNCTION_END(state);
}
