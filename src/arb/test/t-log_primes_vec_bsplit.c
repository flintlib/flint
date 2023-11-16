/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "arb.h"

TEST_FUNCTION_START(arb_log_primes_vec_bsplit, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_ptr v;
        arb_t t;
        slong n, j, prec;
        ulong p, pprev;

        prec = 2 + n_randint(state, 700);
        n = n_randint(state, 30);

        flint_set_num_threads(1 + n_randint(state, 3));

        v = _arb_vec_init(n);
        arb_init(t);

        arb_log_primes_vec_bsplit(v, n, prec);

        pprev = 1;
        for (j = 0; j < n; j++)
        {
            p = n_nextprime(pprev, 1);
            if (p == 2)
                arb_const_log2(t, prec);
            else
                arb_log_ui_from_prev(t, p, t, pprev, prec);
            pprev = p;

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
