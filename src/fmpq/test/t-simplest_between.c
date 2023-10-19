/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_simplest_between, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_t m2, m1, m, l, r, l1, r1;

        fmpq_init(m2);
        fmpq_init(m1);
        fmpq_init(m);
        fmpq_init(l);
        fmpq_init(r);
        fmpq_init(l1);
        fmpq_init(r1);

        fmpq_randtest(l, state, 1 + n_randint(state, 1000));
        fmpq_randtest(r, state, 1 + n_randint(state, 1000));
        fmpz_one(fmpq_numref(m));
        fmpz_randtest_not_zero(fmpq_denref(m), state, 1 + n_randint(state, 1000));
        fmpz_abs(fmpq_denref(m), fmpq_denref(m));
        fmpq_mul(l, l, m);
        fmpq_mul(r, r, m);

        if (i % 2)
        {
            fmpq_randtest(m2, state, 1 + n_randint(state, 1000));
        }
        else
        {
            slong k, bound;
            fmpz * c1;

            bound = 1 + n_randint(state, 100);
            c1 = _fmpz_vec_init(bound);

            fmpz_randtest(c1 + 0, state, 2*FLINT_BITS);
            for (k = 1; k < bound; k++)
            {
                fmpz_randtest_unsigned(c1 + k, state, 3*FLINT_BITS);
                fmpz_add_ui(c1 + k, c1 + k, 1);
            }

            fmpq_set_cfrac(m2, c1, bound);

            _fmpz_vec_clear(c1, bound);
        }

        fmpq_add(l, l, m2);
        fmpq_add(r, r, m2);

        fmpq_simplest_between(m, l, r);

        if (fmpq_cmp(l, r) > 0)
            fmpq_swap(l, r);

        if (!(fmpq_cmp(l, m) <= 0 && fmpq_cmp(m, r) <= 0))
        {
            printf("FAIL\n");
            flint_printf("Check answer is between, i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (   fmpz_cmp(fmpq_denref(m), fmpq_denref(l)) > 0
            || fmpz_cmp(fmpq_denref(m), fmpq_denref(r)) > 0)
        {
            printf("FAIL\n");
            flint_printf("Check denominator, i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (fmpz_cmp_ui(fmpq_denref(m), 1) > 0)
        {
            fmpq_farey_neighbors(l1, r1, m, fmpq_denref(m));
            if (fmpq_cmp(l1, l) >= 0 || fmpq_cmp(r1, r) <= 0)
            {
                printf("FAIL\n");
                flint_printf("Check answer is simplest, i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }
        else if (fmpz_is_one(fmpq_denref(m)))
        {
            fmpz_cdiv_q(fmpq_numref(m1), fmpq_numref(l), fmpq_denref(l));
            if (!fmpz_equal(fmpq_numref(m1), fmpq_numref(m)))
            {
                printf("FAIL\n");
                flint_printf("Check int answer is simplest, i = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            printf("FAIL\n");
            flint_printf("Check return has positive denominator, i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpq_simplest_between(m1, l, m);
        if (!fmpq_equal(m1, m))
        {
            printf("FAIL\n");
            flint_printf("Check m left, i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpq_simplest_between(m1, r, m);
        if (!fmpq_equal(m1, m))
        {
            printf("FAIL\n");
            flint_printf("Check m right, i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(m2);
        fmpq_clear(m1);
        fmpq_clear(m);
        fmpq_clear(l);
        fmpq_clear(r);
        fmpq_clear(l1);
        fmpq_clear(r1);
    }

    TEST_FUNCTION_END(state);
}
