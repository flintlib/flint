/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "bool_mat.h"

TEST_FUNCTION_START(bool_mat_mul, state)
{
    slong iter;

    /* test aliasing */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        bool_mat_t A, B, C, D;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        /* test aliasing with a */
        {
            bool_mat_init(A, m, n);
            bool_mat_init(B, n, n);
            bool_mat_init(C, m, n);
            bool_mat_init(D, m, n);

            bool_mat_randtest(A, state);
            bool_mat_randtest(B, state);
            bool_mat_mul(C, A, B);

            bool_mat_set(D, A);
            bool_mat_mul(D, D, B);

            if (!bool_mat_equal(D, C))
            {
                flint_printf("FAIL (aliasing 1)\n");
                flint_abort();
            }

            bool_mat_clear(A);
            bool_mat_clear(B);
            bool_mat_clear(C);
            bool_mat_clear(D);
        }

        /* test aliasing with b */
        {
            bool_mat_init(A, m, m);
            bool_mat_init(B, m, n);
            bool_mat_init(C, m, n);
            bool_mat_init(D, m, n);

            bool_mat_randtest(A, state);
            bool_mat_randtest(B, state);
            bool_mat_mul(C, A, B);

            bool_mat_set(D, B);
            bool_mat_mul(D, A, D);

            if (!bool_mat_equal(D, C))
            {
                flint_printf("FAIL (aliasing 2)\n");
                flint_abort();
            }

            bool_mat_clear(A);
            bool_mat_clear(B);
            bool_mat_clear(C);
            bool_mat_clear(D);
        }
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, k, l;
        bool_mat_t a, b, c, d, ab, ac, bd, cd, s;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);
        l = n_randint(state, 10);

        bool_mat_init(a, m, n);
        bool_mat_init(b, n, k);
        bool_mat_init(c, n, k);
        bool_mat_init(d, k, l);

        bool_mat_randtest(a, state);
        bool_mat_randtest(b, state);
        bool_mat_randtest(c, state);
        bool_mat_randtest(d, state);

        bool_mat_init(ab, m, k);
        bool_mat_init(ac, m, k);
        bool_mat_init(bd, n, l);
        bool_mat_init(cd, n, l);
        bool_mat_init(s, n, k);

        bool_mat_mul(ab, a, b);
        bool_mat_mul(ac, a, c);
        bool_mat_mul(bd, b, d);
        bool_mat_mul(cd, c, d);
        bool_mat_add(s, b, c);

        /* check associativity of multiplication */
        /* (A*B)*D = A*(B*D) */
        {
            bool_mat_t lhs, rhs;

            bool_mat_init(lhs, m, l);
            bool_mat_init(rhs, m, l);

            bool_mat_mul(lhs, ab, d);
            bool_mat_mul(rhs, a, bd);

            if (!bool_mat_equal(lhs, rhs))
            {
                flint_printf("FAIL (associativity)\n");
                flint_printf("m, n, k, l = %wd, %wd, %wd, %wd\n", m, n, k, l);
                flint_abort();
            }

            bool_mat_clear(lhs);
            bool_mat_clear(rhs);
        }

        /* check left distributivity of multiplication over addition */
        /* A*(B + C) = A*B + A*C */
        {
            bool_mat_t lhs, rhs;

            bool_mat_init(lhs, m, k);
            bool_mat_init(rhs, m, k);

            bool_mat_mul(lhs, a, s);
            bool_mat_add(rhs, ab, ac);

            if (!bool_mat_equal(lhs, rhs))
            {
                flint_printf("FAIL (left distributivity)\n");
                flint_printf("m, n, k, l = %wd, %wd, %wd, %wd\n", m, n, k, l);
                flint_abort();
            }

            bool_mat_clear(lhs);
            bool_mat_clear(rhs);
        }

        /* check right distributivity of multiplication over addition */
        /* (B + C)*D = B*D + C*D */
        {
            bool_mat_t lhs, rhs;

            bool_mat_init(lhs, n, l);
            bool_mat_init(rhs, n, l);

            bool_mat_mul(lhs, s, d);
            bool_mat_add(rhs, bd, cd);

            if (!bool_mat_equal(lhs, rhs))
            {
                flint_printf("FAIL (right distributivity)\n");
                flint_printf("m, n, k, l = %wd, %wd, %wd, %wd\n", m, n, k, l);
                flint_abort();
            }

            bool_mat_clear(lhs);
            bool_mat_clear(rhs);
        }

        /* check left multiplicative identity I*D = D */
        {
            bool_mat_t one, lhs;

            bool_mat_init(one, k, k);
            bool_mat_init(lhs, k, l);

            bool_mat_one(one);
            bool_mat_mul(lhs, one, d);

            if (!bool_mat_equal(lhs, d))
            {
                flint_printf("FAIL (left identity)\n");
                flint_printf("k = %wd, l = %wd\n", k, l);
                flint_abort();
            }

            bool_mat_clear(one);
            bool_mat_clear(lhs);
        }

        /* check right multiplicative identity A*I = A */
        {
            bool_mat_t one, lhs;

            bool_mat_init(one, n, n);
            bool_mat_init(lhs, m, n);

            bool_mat_one(one);
            bool_mat_mul(lhs, a, one);

            if (!bool_mat_equal(lhs, a))
            {
                flint_printf("FAIL (right identity)\n");
                flint_printf("m = %wd, n = %wd\n", m, n);
                flint_abort();
            }

            bool_mat_clear(one);
            bool_mat_clear(lhs);
        }

        bool_mat_clear(a);
        bool_mat_clear(b);
        bool_mat_clear(c);
        bool_mat_clear(d);
        bool_mat_clear(ab);
        bool_mat_clear(ac);
        bool_mat_clear(bd);
        bool_mat_clear(cd);
        bool_mat_clear(s);
    }

    TEST_FUNCTION_END(state);
}
