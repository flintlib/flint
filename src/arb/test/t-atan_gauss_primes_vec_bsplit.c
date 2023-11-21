/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

/* Consecutively, real and imaginary parts of first 64 nonreal
   Gaussian primes. */
static const signed char small_gaussian_primes[] = {
    1, 1, 1, 2, 2, 3, 1, 4, 2, 5, 1, 6, 4, 5, 2, 7, 5, 6, 3, 8, 5, 8, 4, 9,
    1, 10, 3, 10, 7, 8, 4, 11, 7, 10, 6, 11, 2, 13, 9, 10, 7, 12, 1, 14, 2, 15, 8, 13,
    4, 15, 1, 16, 10, 13, 9, 14, 5, 16, 2, 17, 12, 13, 11, 14, 9, 16, 5, 18, 8, 17, 7, 18,
    10, 17, 6, 19, 1, 20, 3, 20, 14, 15, 12, 17, 7, 20, 4, 21, 10, 19, 5, 22, 11, 20, 10, 21,
    14, 19, 13, 20, 1, 24, 8, 23, 5, 24, 17, 18, 16, 19, 4, 25, 13, 22, 6, 25, 12, 23, 1, 26,
    5, 26, 15, 22, 2, 27, 9, 26
};

TEST_FUNCTION_START(arb_atan_gauss_primes_vec_bsplit, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_ptr v;
        arb_t t;
        slong n, j, prec;

        prec = 2 + n_randint(state, 700);
        n = n_randint(state, 30);

        flint_set_num_threads(1 + n_randint(state, 3));

        v = _arb_vec_init(n);
        arb_init(t);

        arb_atan_gauss_primes_vec_bsplit(v, n, prec);

        for (j = 0; j < n; j++)
        {
            arb_set_ui(t, small_gaussian_primes[2 * j + 1]);
            arb_div_ui(t, t, small_gaussian_primes[2 * j], prec);
            arb_atan(t, t, prec);

            if (!arb_overlaps(v + j, t) || arb_rel_accuracy_bits(v + j) < prec - 5)
            {
                flint_printf("FAIL\n\n");
                flint_printf("n = %wu, j = %wd\n", n, j);
                flint_printf("v = "); arb_printd(v + j, 100); flint_printf("\n\n");
                flint_printf("t = "); arb_printd(t, 100); flint_printf("\n\n");
                flint_printf("%wd, %wd\n", prec, arb_rel_accuracy_bits(v + j));
                flint_abort();
            }
        }

        _arb_vec_clear(v, n);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
