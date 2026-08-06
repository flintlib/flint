/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix.h"

static ulong
radix_sub_naive(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    ulong cy, hi, lo, B = LIMB_RADIX(radix);
    slong i;

    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(an >= bn);

    cy = 0;

    for (i = 0; i < an; i++)
    {
        sub_ddmmss(hi, lo, 0, a[i], 0, (i < bn) ? b[i] : 0);
        sub_ddmmss(hi, lo, hi, lo, 0, cy);
        cy = ((hi != 0) || (hi == 0 && lo >= B));
        if (cy)
            lo  += B;
        res[i] = lo;
    }

    return cy;
}

TEST_FUNCTION_START(radix_sub, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c, d;
        slong an, bn;
        ulong cy1, cy2;

        radix_init_randtest(radix, state);

        bn = 1 + n_randint(state, 5);
        an = bn + n_randint(state, 5);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        c = flint_malloc(an * sizeof(ulong));
        d = flint_malloc(an * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn, radix);
        radix_randtest_limbs(c, state, an, radix);
        radix_randtest_limbs(d, state, an, radix);

        cy1 = radix_sub(c, a, an, b, bn, radix);
        cy2 = radix_sub_naive(d, a, an, b, bn, radix);

        if (cy1 != cy2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: sub\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, bn);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", cy1, cy2);
            flint_abort();
        }

        flint_mpn_copyi(c, a, an);
        cy1 = radix_sub(c, c, an, b, bn, radix);

        if (cy1 != cy2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: sub (aliasing a)\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, bn);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", cy1, cy2);
            flint_abort();
        }

        flint_mpn_copyi(c, b, bn);
        cy1 = radix_sub(c, a, an, c, bn, radix);

        if (cy1 != cy2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: sub (aliasing b)\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, bn);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", cy1, cy2);
            flint_abort();
        }

        cy1 = radix_add_naive(c, a, an, b, bn, radix);
        cy1 -= radix_sub_naive(c, c, an, b, bn, radix);

        if (cy1 != 0 || mpn_cmp(a, c, an) != 0)
        {
            flint_printf("FAIL: sub_naive\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, bn);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", cy1, cy2);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
        flint_free(d);
    }

    TEST_FUNCTION_END(state);
}
