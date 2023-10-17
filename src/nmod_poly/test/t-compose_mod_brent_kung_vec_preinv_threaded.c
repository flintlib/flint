/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_compose_mod_brent_kung_vec_preinv_threaded, state)
{
#if FLINT_USES_PTHREAD && (FLINT_USES_TLS || FLINT_REENTRANT)
    int i;
    slong max_threads = 5;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, ainv, b, c;
        mp_limb_t m = n_randtest_prime(state, 0);
        slong j, k, l;
        nmod_poly_struct * pow, * res;

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(ainv, m);

        nmod_poly_randtest(b, state, 1 + n_randint(state, 200));
        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 200));
        l = n_randint(state, 100) + 1;
        k = n_randint(state, l) + 1;

        nmod_poly_rem(b, b, a);
        nmod_poly_reverse(ainv, a, a->length);
        nmod_poly_inv_series(ainv, ainv, a->length);
        pow = (nmod_poly_struct *) flint_malloc((l + k)*sizeof(nmod_poly_struct));
        res = pow + l;

        for (j = 0; j < l; j++)
        {
            nmod_poly_init(pow + j, m);
            nmod_poly_randtest(pow + j, state, n_randint(state, 200) + 1);
            nmod_poly_rem(pow + j, pow + j, a);
        }

        for (j = 0; j < k; j++)
            nmod_poly_init(res + j, m);

        nmod_poly_compose_mod_brent_kung_vec_preinv_threaded(res, pow, l, k,
                b, a, ainv);

        for (j = 0; j < k; j++)
        {
            nmod_poly_compose_mod(c, pow + j, b, a);
            if (!nmod_poly_equal(res + j, c))
            {
                flint_printf("FAIL (composition):\n");
                flint_printf("a:\n"); nmod_poly_print(a); flint_printf("\n");
                flint_printf("res:\n"); nmod_poly_print(res + j); flint_printf("\n");
                flint_printf("pow:\n"); nmod_poly_print(pow + j); flint_printf("\n");
                flint_printf("b:\n"); nmod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); nmod_poly_print(c); flint_printf("\n");
                flint_printf("j: %wd\n", j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(ainv);
        nmod_poly_clear(b);
        nmod_poly_clear(c);

        for (j = 0; j < l; j++)
            nmod_poly_clear(pow + j);

        for (j = 0; j < k; j++)
            nmod_poly_clear(res + j);

        flint_free(pow);
    }

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
