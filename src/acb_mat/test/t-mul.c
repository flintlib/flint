/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
_acb_mat_init_randtest(acb_mat_t mat, slong r, slong c, flint_rand_t state)
{
    acb_mat_init(mat, r, c);
    acb_mat_randtest(mat, state, 2 + n_randint(state, 200), 10);
}

void
_acb_mat_nprintd(const char * name, acb_mat_t mat)
{
    flint_printf("%s = ", name);
    acb_mat_printd(mat, 15);
    flint_printf("\n\n");
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong m, n, k, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C;
        acb_mat_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, k);
        fmpq_mat_init(C, m, k);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, k);
        acb_mat_init(c, m, k);
        acb_mat_init(d, m, k);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_randtest(B, state, qbits2);
        fmpq_mat_mul(C, A, B);

        acb_mat_set_fmpq_mat(a, A, rbits1);
        acb_mat_set_fmpq_mat(b, B, rbits2);
        acb_mat_mul(c, a, b, rbits3);

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
            acb_mat_mul(d, d, b, rbits3);
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
            acb_mat_mul(d, a, d, rbits3);
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

    /* check algebraic properties like associativity and distributivity */
    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        slong m, n, k, l;
        slong rbits;
        acb_mat_t a, b, c, d, ab, ac, bd, cd, s;

        rbits = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);
        l = n_randint(state, 10);

        _acb_mat_init_randtest(a, m, n, state);
        _acb_mat_init_randtest(b, n, k, state);
        _acb_mat_init_randtest(c, n, k, state);
        _acb_mat_init_randtest(d, k, l, state);

        acb_mat_init(ab, m, k);
        acb_mat_init(ac, m, k);
        acb_mat_init(bd, n, l);
        acb_mat_init(cd, n, l);
        acb_mat_init(s, n, k);

        acb_mat_mul(ab, a, b, rbits);
        acb_mat_mul(ac, a, c, rbits);
        acb_mat_mul(bd, b, d, rbits);
        acb_mat_mul(cd, c, d, rbits);
        acb_mat_add(s, b, c, rbits);

        /* check associativity of multiplication */
        /* (A*B)*D = A*(B*D) */
        {
            acb_mat_t lhs, rhs;

            acb_mat_init(lhs, m, l);
            acb_mat_init(rhs, m, l);

            acb_mat_mul(lhs, ab, d, rbits);
            acb_mat_mul(rhs, a, bd, rbits);

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

            acb_mat_init(lhs, m, k);
            acb_mat_init(rhs, m, k);

            acb_mat_mul(lhs, a, s, rbits);
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

            acb_mat_init(lhs, n, l);
            acb_mat_init(rhs, n, l);

            acb_mat_mul(lhs, s, d, rbits);
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

            acb_mat_init(one, k, k);
            acb_mat_init(lhs, k, l);

            acb_mat_one(one);
            acb_mat_mul(lhs, one, d, rbits);

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

            acb_mat_init(one, n, n);
            acb_mat_init(lhs, m, n);

            acb_mat_one(one);
            acb_mat_mul(lhs, a, one, rbits);

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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
