/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("mullow_classical....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpq_poly */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong qbits1, qbits2, rbits1, rbits2, rbits3, trunc;
        fmpq_poly_t A, B, C;
        acb_poly_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);
        trunc = n_randint(state, 10);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(C);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        fmpq_poly_randtest(A, state, 1 + n_randint(state, 10), qbits1);
        fmpq_poly_randtest(B, state, 1 + n_randint(state, 10), qbits2);
        fmpq_poly_mullow(C, A, B, trunc);

        acb_poly_set_fmpq_poly(a, A, rbits1);
        acb_poly_set_fmpq_poly(b, B, rbits2);
        acb_poly_mullow_classical(c, a, b, trunc, rbits3);

        if (!acb_poly_contains_fmpq_poly(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits3 = %wd\n", rbits3);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_poly_print(C); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_mullow_classical(d, d, b, trunc, rbits3);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_abort();
        }

        acb_poly_set(d, b);
        acb_poly_mullow_classical(d, a, d, trunc, rbits3);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
            flint_abort();
        }

        /* test squaring */
        acb_poly_set(b, a);
        acb_poly_mullow_classical(c, a, b, trunc, rbits3);
        acb_poly_mullow_classical(d, a, a, trunc, rbits3);
        if (!acb_poly_overlaps(c, d))  /* not guaranteed to be identical */
        {
            flint_printf("FAIL (squaring)\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_mullow_classical(a, a, a, trunc, rbits3);
        if (!acb_poly_equal(d, a))
        {
            flint_printf("FAIL (aliasing, squaring)\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    /* check a*(b+c) = a*b+a*c */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong bits, trunc;
        acb_poly_t a, b, c, bc, abc, ab, ac, abac;

        bits = 2 + n_randint(state, 200);
        trunc = n_randint(state, 10);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(bc);
        acb_poly_init(abc);
        acb_poly_init(ab);
        acb_poly_init(ac);
        acb_poly_init(abac);

        acb_poly_randtest(a, state, 1 + n_randint(state, 10), bits, 5);
        acb_poly_randtest(b, state, 1 + n_randint(state, 10), bits, 5);
        acb_poly_randtest(c, state, 1 + n_randint(state, 10), bits, 5);

        acb_poly_add(bc, b, c, bits);
        acb_poly_mullow_classical(abc, a, bc, trunc, bits);

        acb_poly_mullow_classical(ab, a, b, trunc, bits);
        acb_poly_mullow_classical(ac, a, c, trunc, bits);
        acb_poly_add(abac, ab, ac, bits);

        if (!acb_poly_overlaps(abc, abac))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("abc = "); acb_poly_printd(abc, 15); flint_printf("\n\n");
            flint_printf("abac = "); acb_poly_printd(abac, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(bc);
        acb_poly_clear(abc);
        acb_poly_clear(ab);
        acb_poly_clear(ac);
        acb_poly_clear(abac);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

