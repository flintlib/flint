/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mulmod_precond, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res1, res2, t, d, dinv;
        nmod_poly_mulmod_precond_t precond;
        int method;

        ulong p = n_randtest_prime(state, 0);

        nmod_t mod;
        nmod_init(&mod, p);

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(d, p);
        nmod_poly_init(res1, p);
        nmod_poly_init(res2, p);
        nmod_poly_init(t, p);
        nmod_poly_init(dinv, p);

        switch (n_randint(state, 3))
        {
            case 0: method = NMOD_POLY_MULMOD_PRECOND_NONE; break;
            case 1: method = NMOD_POLY_MULMOD_PRECOND_SHOUP; break;
            case 2: method = NMOD_POLY_MULMOD_PRECOND_MATRIX; break;
            default: flint_abort();
        }

        do {
            nmod_poly_randtest(d, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(d));

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));

        if (a->length >= d->length)
          nmod_poly_rem(a, a, d);
        if (b->length >= d->length)
          nmod_poly_rem(b, b, d);

        nmod_poly_randtest(res1, state, n_randint(state, 50));

        nmod_poly_reverse(dinv, d, d->length);
        nmod_poly_inv_series(dinv, dinv, d->length);

        if (n_randint(state, 2))
            nmod_poly_mulmod_precond_init_method(precond, a, d, dinv, method);
        else
            nmod_poly_mulmod_precond_init_num(precond, a, d, dinv, n_randint(state, 100));

        nmod_poly_mulmod_precond(res1, precond, b);

        nmod_poly_mul(res2, a, b);
        nmod_poly_divrem(t, res2, res2, d);

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("p = %wu, method = %d\n\n", p, precond->method);
            flint_printf("a = %{nmod_poly}\n\n", a);
            flint_printf("b = %{nmod_poly}\n\n", b);
            flint_printf("d = %{nmod_poly}\n\n", d);
            flint_printf("res1 = %{nmod_poly}\n\n", res1);
            flint_printf("res2 = %{nmod_poly}\n\n", res2);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mulmod_precond_clear(precond);

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(d);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
        nmod_poly_clear(dinv);
    }

    TEST_FUNCTION_END(state);
}
