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
#include "fmpz.h"

static ulong
radix_divrem_1_naive(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
{
    ulong B = LIMB_RADIX(radix);
    slong i;
    ulong r;

    fmpz_t T, D, Q, R;

    fmpz_init(T);
    fmpz_init(D);
    fmpz_init(Q);
    fmpz_init(R);

    for (i = 0; i < an; i++)
    {
        fmpz_ui_pow_ui(Q, B, i);
        fmpz_addmul_ui(T, Q, a[i]);
    }

    fmpz_set_ui(D, d);

    fmpz_fdiv_qr(Q, R, T, D);

    for (i = 0; i < an; i++)
    {
        res[i] = fmpz_fdiv_ui(Q, B);
        fmpz_fdiv_q_ui(Q, Q, B);
    }

    r = fmpz_get_ui(R);

    fmpz_clear(T);
    fmpz_clear(D);
    fmpz_clear(Q);
    fmpz_clear(R);

    return r;
}

TEST_FUNCTION_START(radix_divrem_1, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, c, d;
        ulong b;
        slong an;
        ulong r1, r2;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 5);

        a = flint_malloc(an * sizeof(ulong));
        c = flint_malloc(an * sizeof(ulong));
        d = flint_malloc(an * sizeof(ulong));

        b = n_randint(state, LIMB_RADIX(radix));
        b = FLINT_MAX(b, 1);

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(c, state, an, radix);
        radix_randtest_limbs(d, state, an, radix);

        r1 = radix_divrem_1(c, a, an, b, radix);
        r2 = radix_divrem_1_naive(d, a, an, b, radix);

        if (r1 != r2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: divrem_1\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%wu\n", b);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", r1, r2);
            flint_abort();
        }

        flint_mpn_copyi(c, a, an);
        r1 = radix_divrem_1(c, c, an, b, radix);

        if (r1 != r2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: divrem_1 (aliasing)\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%wu\n", b);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", r1, r2);
            flint_abort();
        }

        r1 = radix_divrem_two(c, a, an, radix);
        r2 = radix_divrem_1_naive(d, a, an, 2, radix);

        if (r1 != r2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: divrem_two\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", r1, r2);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(c);
        flint_free(d);
    }

    TEST_FUNCTION_END(state);
}
