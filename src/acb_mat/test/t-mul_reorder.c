/*
    Copyright (C) 2012, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

/* Defined in t-mul.c and t-mul_reorder.c */
#ifndef _acb_mat_init_randtest
#define _acb_mat_init_randtest _acb_mat_init_randtest
void _acb_mat_init_randtest(acb_mat_t mat, slong r, slong c, flint_rand_t state)
{
    acb_mat_init(mat, r, c);
    acb_mat_randtest(mat, state, 2 + n_randint(state, 200), 10);
}
#endif

/* Defined in t-mul.c and t-mul_reorder.c */
#ifndef _acb_mat_nprintd
#define _acb_mat_nprintd _acb_mat_nprintd
void
_acb_mat_nprintd(const char * name, acb_mat_t mat)
{
    flint_printf("%s = ", name);
    acb_mat_printd(mat, 15);
    flint_printf("\n\n");
}
#endif

TEST_FUNCTION_START(acb_mat_mul_reorder, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, k, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C;
        acb_mat_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = n_randint(state, 8);
        n = n_randint(state, 8);
        k = n_randint(state, 8);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, k);
        fmpq_mat_init(C, m, k);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, k);
        _acb_mat_init_randtest(c, m, k, state);
        _acb_mat_init_randtest(d, m, k, state);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_randtest(B, state, qbits2);
        fmpq_mat_mul(C, A, B);

        acb_mat_set_fmpq_mat(a, A, rbits1);
        acb_mat_set_fmpq_mat(b, B, rbits2);
        acb_mat_mul_reorder(c, a, b, rbits3);

        if (!acb_mat_contains_fmpq_mat(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd, k = %wd, bits3 = %wd\n", m, n, k, rbits3);

            flint_printf("A = "); fmpq_mat_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_mat_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_mat_print(C); flint_printf("\n\n");

            flint_printf("a = "); acb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_mat_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_mat_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* test aliasing with a */
        if (acb_mat_nrows(a) == acb_mat_nrows(c) &&
            acb_mat_ncols(a) == acb_mat_ncols(c))
        {
            acb_mat_set(d, a);
            acb_mat_mul_reorder(d, d, b, rbits3);
            if (!acb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 1)\n\n");
                flint_abort();
            }
        }

        /* test aliasing with b */
        if (acb_mat_nrows(b) == acb_mat_nrows(c) &&
            acb_mat_ncols(b) == acb_mat_ncols(c))
        {
            acb_mat_set(d, b);
            acb_mat_mul_reorder(d, a, d, rbits3);
            if (!acb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 2)\n\n");
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);
    }

    /* general aliasing test */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        slong rbits;
        acb_mat_t a, b, c, d;

        rbits = 2 + n_randint(state, 200);

        m = n_randint(state, 8);
        n = n_randint(state, 8);

        _acb_mat_init_randtest(a, m, n, state);
        _acb_mat_init_randtest(b, n, n, state);
        _acb_mat_init_randtest(c, m, n, state);
        _acb_mat_init_randtest(d, m, n, state);

        acb_mat_mul_reorder(c, a, b, rbits);
        acb_mat_set(d, a);
        acb_mat_mul_reorder(d, d, b, rbits);

        if (!acb_mat_equal(c, d))
        {
            flint_printf("FAIL (aliasing 3)\n\n");
            _acb_mat_nprintd("a", a);
            _acb_mat_nprintd("b", b);
            _acb_mat_nprintd("c", d);
            _acb_mat_nprintd("d", d);
            flint_abort();
        }

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);

        _acb_mat_init_randtest(a, m, m, state);
        _acb_mat_init_randtest(b, m, n, state);
        _acb_mat_init_randtest(c, m, n, state);
        _acb_mat_init_randtest(d, m, n, state);

        acb_mat_mul_reorder(c, a, b, rbits);
        acb_mat_set(d, b);
        acb_mat_mul_reorder(d, a, d, rbits);

        if (!acb_mat_equal(c, d))
        {
            flint_printf("FAIL (aliasing 4)\n\n");
            _acb_mat_nprintd("a", a);
            _acb_mat_nprintd("b", b);
            _acb_mat_nprintd("c", d);
            _acb_mat_nprintd("d", d);
            flint_abort();
        }

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);
    }

    /* check algebraic properties like associativity and distributivity */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, k, l;
        slong rbits;
        acb_mat_t a, b, c, d, ab, ac, bd, cd, s;

        rbits = 2 + n_randint(state, 200);

        m = n_randint(state, 8);
        n = n_randint(state, 8);
        k = n_randint(state, 8);
        l = n_randint(state, 8);

        _acb_mat_init_randtest(a, m, n, state);
        _acb_mat_init_randtest(b, n, k, state);
        _acb_mat_init_randtest(c, n, k, state);
        _acb_mat_init_randtest(d, k, l, state);

        _acb_mat_init_randtest(ab, m, k, state);
        _acb_mat_init_randtest(ac, m, k, state);
        _acb_mat_init_randtest(bd, n, l, state);
        _acb_mat_init_randtest(cd, n, l, state);
        _acb_mat_init_randtest(s, n, k, state);

        acb_mat_mul_reorder(ab, a, b, rbits);
        acb_mat_mul_reorder(ac, a, c, rbits);
        acb_mat_mul_reorder(bd, b, d, rbits);
        acb_mat_mul_reorder(cd, c, d, rbits);
        acb_mat_add(s, b, c, rbits);

        /* check associativity of multiplication */
        /* (A*B)*D = A*(B*D) */
        {
            acb_mat_t lhs, rhs;

            _acb_mat_init_randtest(lhs, m, l, state);
            _acb_mat_init_randtest(rhs, m, l, state);

            acb_mat_mul_reorder(lhs, ab, d, rbits);
            acb_mat_mul_reorder(rhs, a, bd, rbits);

            if (!acb_mat_overlaps(lhs, rhs))
            {
                flint_printf("FAIL\n\n");
                flint_printf("m, n, k, l = %wd, %wd, %wd, %wd\n", m, n, k, l);
                flint_printf("rbits = %wd\n", rbits);

                _acb_mat_nprintd("a", a);
                _acb_mat_nprintd("b", b);
                _acb_mat_nprintd("d", d);
                _acb_mat_nprintd("(a*b)*d", lhs);
                _acb_mat_nprintd("a*(b*d)", rhs);

                flint_abort();
            }

            acb_mat_clear(lhs);
            acb_mat_clear(rhs);
        }

        /* check left distributivity of multiplication over addition */
        /* A*(B + C) = A*B + A*C */
        {
            acb_mat_t lhs, rhs;

            _acb_mat_init_randtest(lhs, m, k, state);
            _acb_mat_init_randtest(rhs, m, k, state);

            acb_mat_mul_reorder(lhs, a, s, rbits);
            acb_mat_add(rhs, ab, ac, rbits);

            if (!acb_mat_overlaps(lhs, rhs))
            {
                flint_printf("FAIL\n\n");
                flint_printf("m, n, k, l = %wd, %wd, %wd, %wd\n", m, n, k, l);
                flint_printf("rbits = %wd\n", rbits);

                _acb_mat_nprintd("a", a);
                _acb_mat_nprintd("b", b);
                _acb_mat_nprintd("c", c);
                _acb_mat_nprintd("a*(b + c)", lhs);
                _acb_mat_nprintd("a*b + b*c", rhs);

                flint_abort();
            }

            acb_mat_clear(lhs);
            acb_mat_clear(rhs);
        }

        /* check right distributivity of multiplication over addition */
        /* (B + C)*D = B*D + C*D */
        {
            acb_mat_t lhs, rhs;

            _acb_mat_init_randtest(lhs, n, l, state);
            _acb_mat_init_randtest(rhs, n, l, state);

            acb_mat_mul_reorder(lhs, s, d, rbits);
            acb_mat_add(rhs, bd, cd, rbits);

            if (!acb_mat_overlaps(lhs, rhs))
            {
                flint_printf("FAIL\n\n");
                flint_printf("m, n, k, l = %wd, %wd, %wd, %wd\n", m, n, k, l);
                flint_printf("rbits = %wd\n", rbits);

                _acb_mat_nprintd("b", b);
                _acb_mat_nprintd("c", c);
                _acb_mat_nprintd("d", d);
                _acb_mat_nprintd("(b + c)*d", lhs);
                _acb_mat_nprintd("b*d + c*d", rhs);

                flint_abort();
            }

            acb_mat_clear(lhs);
            acb_mat_clear(rhs);
        }

        /* check left multiplicative identity I*D = D */
        {
            acb_mat_t one, lhs;

            _acb_mat_init_randtest(one, k, k, state);
            _acb_mat_init_randtest(lhs, k, l, state);

            acb_mat_one(one);
            acb_mat_mul_reorder(lhs, one, d, rbits);

            if (!acb_mat_contains(lhs, d))
            {
                flint_printf("FAIL\n\n");
                flint_printf("k = %wd, l = %wd\n", k, l);
                flint_printf("rbits = %wd\n", rbits);

                _acb_mat_nprintd("identity * d", lhs);
                _acb_mat_nprintd("d", d);

                flint_abort();
            }

            acb_mat_clear(one);
            acb_mat_clear(lhs);
        }

        /* check right multiplicative identity A*I = A */
        {
            acb_mat_t one, lhs;

            _acb_mat_init_randtest(one, n, n, state);
            _acb_mat_init_randtest(lhs, m, n, state);

            acb_mat_one(one);
            acb_mat_mul_reorder(lhs, a, one, rbits);

            if (!acb_mat_contains(lhs, a))
            {
                flint_printf("FAIL\n\n");
                flint_printf("m = %wd, n = %wd\n", m, n);
                flint_printf("rbits = %wd\n", rbits);

                _acb_mat_nprintd("a * identity", lhs);
                _acb_mat_nprintd("a", a);

                flint_abort();
            }

            acb_mat_clear(one);
            acb_mat_clear(lhs);
        }

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);
        acb_mat_clear(ab);
        acb_mat_clear(ac);
        acb_mat_clear(bd);
        acb_mat_clear(cd);
        acb_mat_clear(s);
    }

    TEST_FUNCTION_END(state);
}
