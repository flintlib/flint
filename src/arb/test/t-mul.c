/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

#ifndef mag_close
#define mag_close mag_close
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
#endif

#ifndef arb_equal_mid_close_mag
#define arb_equal_mid_close_mag arb_equal_mid_close_mag
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
#endif

void
arb_mul_naive(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arf_t zm_exact, zm_rounded, zr, t, u;

    arf_init(zm_exact);
    arf_init(zm_rounded);
    arf_init(zr);
    arf_init(t);
    arf_init(u);

    arf_mul(zm_exact, arb_midref(x), arb_midref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_set_round(zm_rounded, zm_exact, prec, ARB_RND);

    /* rounding error */
    if (arf_equal(zm_exact, zm_rounded))
    {
        arf_zero(zr);
    }
    else
    {
        fmpz_t e;
        fmpz_init(e);

        /* more accurate, but not what we are testing
        arf_sub(zr, zm_exact, zm_rounded, MAG_BITS, ARF_RND_UP);
        arf_abs(zr, zr); */

        fmpz_sub_ui(e, ARF_EXPREF(zm_rounded), prec);
        arf_one(zr);
        arf_mul_2exp_fmpz(zr, zr, e);
        fmpz_clear(e);
    }

    /* propagated error - note that (signed infinity) * nonzero should be a signed
       infinity, meaning the error should *not* propagate */
    if (!arb_is_exact(x) && !(
            arb_is_nonzero(x) && arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y))))
    {
        arf_set_mag(t, arb_radref(x));
        arf_abs(u, arb_midref(y));
        arf_addmul(zr, t, u, MAG_BITS, ARF_RND_UP);
    }

    if (!arb_is_exact(y) && !(
            arb_is_nonzero(y) && arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x))))
    {
        arf_set_mag(t, arb_radref(y));
        arf_abs(u, arb_midref(x));
        arf_addmul(zr, t, u, MAG_BITS, ARF_RND_UP);
    }

    if (!arb_is_exact(x) && !arb_is_exact(y))
    {
        arf_set_mag(t, arb_radref(x));
        arf_set_mag(u, arb_radref(y));
        arf_addmul(zr, t, u, MAG_BITS, ARF_RND_UP);
    }

    arf_set(arb_midref(z), zm_rounded);
    arf_get_mag(arb_radref(z), zr);

    arf_clear(zm_exact);
    arf_clear(zm_rounded);
    arf_clear(zr);
    arf_clear(t);
    arf_clear(u);
}

TEST_FUNCTION_START(arb_mul, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
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

        arb_mul(c, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

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
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
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

        arb_mul(a, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

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
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
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

        arb_mul(b, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

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
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
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

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                arb_mul(z, x, y, prec);
                arb_mul_naive(v, x, y, prec);

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
                arb_mul(z, x, y, prec);
                arb_mul(v, x, x, prec);

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
                arb_mul(v, x, x, prec);
                arb_mul(x, x, x, prec);

                if (!arb_equal(v, x))
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 3:
                arb_mul(v, x, y, prec);
                arb_mul(x, x, y, prec);

                if (!arb_equal(v, x))
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            default:
                arb_mul(v, x, y, prec);
                arb_mul(x, y, x, prec);

                if (!arf_equal(arb_midref(x), arb_midref(v))
                    || !mag_close(arb_radref(x), arb_radref(v)))
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

    TEST_FUNCTION_END(state);
}
