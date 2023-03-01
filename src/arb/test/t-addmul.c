/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
mag_close(const mag_t am, const mag_t bm)
{
    arf_t t, a, b;
    int res1, res2;

    arf_init(t);
    arf_init(a);
    arf_init(b);

    arf_set_mag(a, am);
    arf_set_mag(b, bm);

    arf_mul_ui(t, b, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res1 = arf_cmp(a, t) <= 0;

    arf_mul_ui(t, a, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res2 = arf_cmp(b, t) <= 0;

    arf_clear(t);
    arf_clear(a);
    arf_clear(b);

    return res1 && res2;
}

int
arb_equal_mid_close_mag(const arb_t a, const arb_t b)
{
    return arf_equal(arb_midref(a), arb_midref(b)) &&
        (mag_close(arb_radref(a), arb_radref(b)) ||
            /* If a's and b's centres are infinite but their radii are finite, the radii don't need
               to be close: they represent signed infinity regardless. If their centres are NaN,
               then we should ignore their radii. */
            (arf_is_inf(arb_midref(a)) && arf_is_inf(arb_midref(b)) &&
                mag_is_finite(arb_radref(a)) && mag_is_finite(arb_radref(b))) ||
            (arf_is_nan(arb_midref(a)) && arf_is_nan(arb_midref(b))));
}

void
arb_addmul_naive(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_mul(t, x, y, ARF_PREC_EXACT);
    arb_add(z, z, t, prec);
    arb_clear(t);
}

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("addmul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        arb_randtest(c, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(z, state, c, 1 + n_randint(state, 200));

        arb_addmul(c, a, b, 2 + n_randint(state, 200));
        fmpq_addmul(z, x, y);

        if (!arb_contains_fmpq(c, z))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        fmpq_set(z, x);

        arb_addmul(a, a, b, 2 + n_randint(state, 200));
        fmpq_addmul(z, x, y);

        if (!arb_contains_fmpq(a, z))
        {
            flint_printf("FAIL: aliasing (c, a)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        fmpq_set(z, y);

        arb_addmul(b, a, b, 2 + n_randint(state, 200));
        fmpq_addmul(z, x, y);

        if (!arb_contains_fmpq(b, z))
        {
            flint_printf("FAIL: aliasing (c, b)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* main test */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y, z, v;
        slong prec;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);
            arb_randtest_special(y, state, n_randint(state,2) ? 2000 : 200, 200);
            arb_randtest_special(z, state, n_randint(state,2) ? 2000 : 200, 200);

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                arb_set(v, z);
                arb_addmul(z, x, y, prec);
                arb_addmul_naive(v, x, y, prec);

                if (!arb_equal_mid_close_mag(z, v))
                {
                    flint_printf("FAIL!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 1:
                arb_set(y, x);
                arb_set(z, v);
                arb_addmul(z, x, y, prec);
                arb_addmul(v, x, x, prec);

                if (!arb_equal_mid_close_mag(z, v))
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 2:
                arb_set(v, x);
                arb_addmul(v, x, x, prec);
                arb_addmul(x, x, x, prec);

                if (!arb_equal_mid_close_mag(x, v))
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 3:
                arb_set(v, x);
                arb_addmul(v, x, y, prec);
                arb_addmul(x, x, y, prec);

                if (!arb_equal_mid_close_mag(x, v))
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            default:
                arb_set(v, x);
                arb_addmul(v, x, y, prec);
                arb_addmul(x, y, x, prec);

                if (!arb_equal_mid_close_mag(x, v))
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;
            }
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
