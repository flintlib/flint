/*
    Copyright (C) 2022 Erik Postma

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

#define ASSERT(cond) if (!(cond)) { flint_printf("FAIL: %d\n", __LINE__); flint_abort(); }

void naive_nonnegative_abs(arb_t y, const arb_t x)
{
    if(y == x)
    {
        /* We don't want to worry about aliasing */
        arb_t z;
        arb_init(z);
        naive_nonnegative_abs(z, x);
        arb_swap(y, z);
        arb_clear(z);
        return;
    }

    if(! arb_is_finite(x))
    {
        if(arf_cmp_si(arb_midref(x), 0) < 0) {
            /* Negative infinity, or "whole real line" represented with a negative centre. We can't
               represent the correct half interval we would like to return for the latter case; in
               order to match arb_abs and arb_nonnegative_abs, we flip the sign for the centre in
               both cases. */
            arb_neg(y, x);
        }
        else
        {
            /* Positive infinity, whole real line represented with a negative centre, or NaN. */
            arb_set(y, x);
        }

        return;
    }

    if(arf_cmp_si(arb_midref(x), 0) > 0)
    {
        if(arf_cmpabs_mag(arb_midref(x), arb_radref(x)) < 0)
        {
            arf_set_mag(arb_midref(y), arb_radref(x));
            arf_add(arb_midref(y), arb_midref(y), arb_midref(x), MAG_BITS+10, ARF_RND_UP);
            arf_mul_2exp_si(arb_midref(y), arb_midref(y), -1);
            arf_get_mag(arb_radref(y), arb_midref(y));
            arf_set_mag(arb_midref(y), arb_radref(y));
        }
        else
        {
            arb_set(y, x);
        }
    }
    else if(arf_cmp_si(arb_midref(x), 0) < 0)
    {
        if(arf_cmpabs_mag(arb_midref(x), arb_radref(x)) < 0)
        {
            arf_set_mag(arb_midref(y), arb_radref(x));
            arf_sub(arb_midref(y), arb_midref(y), arb_midref(x), MAG_BITS+10, ARF_RND_UP);
            arf_mul_2exp_si(arb_midref(y), arb_midref(y), -1);
            arf_get_mag(arb_radref(y), arb_midref(y));
            arf_set_mag(arb_midref(y), arb_radref(y));
        }
        else
        {
            arb_neg(y, x);
        }
    }
    else
    {
        ASSERT(arf_equal_si(arb_midref(x), 0));
        mag_mul_2exp_si(arb_radref(y), arb_radref(x), -1);
        arf_set_mag(arb_midref(y), arb_radref(y));
    }
}

/* Let u equal 1 ulp of y's radius. Are x's centre and radius within u of y's centre and radius,
   respectively? */
int nearly_equal(const arb_t x, const arb_t y)
{
    arf_t s, t;
    mag_t ulp;
    int res;

    arf_init(t);
    arf_init(s);
    mag_init(ulp);

    arf_set_mag(t, arb_radref(y));
    if(mag_is_special(arb_radref(y)))
    {
        mag_zero(ulp);
    }
    else
    {
        arf_mag_set_ulp(ulp, t, MAG_BITS);
    }
    arf_set_mag(s, arb_radref(x));
    arf_sub(t, t, s, MAG_BITS, ARF_RND_UP);
    res = mag_equal(arb_radref(x), arb_radref(y)) || (arf_cmpabs_mag(t, ulp) <= 0);

    arf_sub(t, arb_midref(x), arb_midref(y), MAG_BITS+2, ARF_RND_UP);
    res = res && (arf_equal(arb_midref(x), arb_midref(y)) || (arf_cmpabs_mag(t, ulp) <= 0));

    mag_clear(ulp);
    arf_clear(s);
    arf_clear(t);

    return res;
}

TEST_FUNCTION_START(arb_nonnegative_abs, state)
{
    slong iter, wide;

    for (wide = 0; wide < 2; wide++)
    {
        for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
        {
            arb_t a, b, c;

            arb_init(a);
            arb_init(b);
            arb_init(c);

            if(wide)
            {
                arb_randtest_wide(a, state, 1 + n_randint(state, 2000), 100);
            }
            else
            {
                arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
            }

            arb_nonnegative_abs(b, a);
            naive_nonnegative_abs(c, a);

            if(! nearly_equal(b, c))
            {
                flint_printf("FAIL: nonnegative_abs\n\n");
                flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 100); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
                flint_abort();
            }

            arb_set(b, a);
            arb_nonnegative_abs(b, b);
            if(! nearly_equal(b, c))
            {
                flint_printf("FAIL: aliasing (nonnegative_abs)\n\n");
                flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 100); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
                flint_abort();
            }

            arb_clear(c);
            arb_clear(b);
            arb_clear(a);
        }
    }

    TEST_FUNCTION_END(state);
}
