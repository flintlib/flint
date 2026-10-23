/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_get_cfrac, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_t x, r;
        fmpz *c1, *c2;
        slong k, n1, n2, bound;

        fmpq_init(x);
        fmpq_init(r);

        if (i % 2)
        {
            fmpq_randtest(x, state, 1 + n_randint(state, 1000));
            bound = fmpq_cfrac_bound(x);
            c1 = _fmpz_vec_init(bound);
        }
        else
        {
            bound = 1 + n_randint(state, 100);
            c1 = _fmpz_vec_init(bound);

            fmpz_randtest(c1 + 0, state, 2*FLINT_BITS);
            for (k = 1; k < bound; k++)
            {
                fmpz_randtest_unsigned(c1 + k, state, 3*FLINT_BITS);
                fmpz_add_ui(c1 + k, c1 + k, 1);
            }

            fmpq_set_cfrac(x, c1, bound);
        }

        c2 = _fmpz_vec_init(bound);

        n1 = fmpq_get_cfrac(c1, r, x, bound);

        if (!fmpq_is_zero(r))
        {
            flint_printf("FAIL: expected zero remainder\n");
            fflush(stdout);
            flint_abort();
        }

        /* Test chaining */
        n2 = 0;
        while (1)
        {
            n2 += fmpq_get_cfrac(c2 + n2, x, x, n_randint(state, 10));
            if (fmpq_is_zero(x))
                break;
            fmpq_inv(x, x);
        }

        if (n1 != n2)
        {
            flint_printf("FAIL: i = %wd, n1 = %wd, n2 = %wd\n", i, n1, n2);
            fflush(stdout);
            flint_abort();
        }

        if (!_fmpz_vec_equal(c1, c2, n1))
        {
            flint_printf("FAIL: i = %wd, vectors not equal\n", i);
            _fmpz_vec_print(c1, n1); flint_printf("\n");
            _fmpz_vec_print(c2, n2); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(c1, bound);
        _fmpz_vec_clear(c2, bound);
        fmpq_clear(x);
        fmpq_clear(r);
    }

    /* No terms: the remainder is 1/x in canonical form */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_t x, r, t;
        fmpz c[1];
        slong n;

        fmpq_init(x);
        fmpq_init(r);
        fmpq_init(t);
        fmpz_init(c);

        if (i == 0)
            fmpq_set_si(x, -3, 2);
        else
            fmpq_randtest_not_zero(x, state, 1 + n_randint(state, 200));

        n = fmpq_get_cfrac(c, r, x, -(slong) n_randint(state, 2));
        fmpq_inv(t, x);

        if (n != 0 || !fmpq_is_canonical(r) || !fmpq_equal(r, t))
        {
            flint_printf("FAIL: no terms\n");
            flint_printf("x = "); fmpq_print(x);
            flint_printf(", r = "); fmpz_print(fmpq_numref(r));
            flint_printf("/"); fmpz_print(fmpq_denref(r));
            flint_printf(", n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(r);
        fmpq_clear(t);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
