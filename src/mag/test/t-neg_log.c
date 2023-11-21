/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arf.h"
#include "mag.h"

/* Defined in t-log.c and t-neg_log.c */
#ifndef arf_log
#define arf_log arf_log
void
arf_log(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
{
    _arf_call_mpfr_func(y, NULL, (int (*)(void)) mpfr_log, x, NULL, prec, rnd);
}
#endif

TEST_FUNCTION_START(mag_neg_log, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, z2;
        fmpz_t n;
        mag_t xb, yb;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(z2);

        mag_init(xb);
        mag_init(yb);

        fmpz_init(n);

        mag_randtest_special(xb, state, 25);
        mag_randtest_special(yb, state, 25);

        if (n_randint(state, 2))
        {
            mag_neg_log(yb, xb);

            arf_set_mag(x, xb);
            if (mag_cmp_2exp_si(xb, 0) >= 0)
                arf_zero(z);
            else
            {
                arf_log(z, x, MAG_BITS, ARF_RND_UP);
                arf_neg(z, z);
            }
        }
        else
        {
            arf_set_mag(x, xb);

            fmpz_randtest(n, state, 100);
            mag_mul_2exp_fmpz(xb, xb, n);
            mag_neg_log(yb, xb);

            if (mag_cmp_2exp_si(xb, 0) >= 0)
            {
                arf_zero(z);
            }
            else
            {
                arf_log(z, x, 2 * MAG_BITS, ARF_RND_UP);
                arf_neg(z, z);
                arf_set_ui(z2, 2);
                arf_log(z2, z2, 2 * MAG_BITS, ARF_RND_UP);
                arf_submul_fmpz(z, z2, n, 2 * MAG_BITS, ARF_RND_UP);
            }
        }

        arf_set_mag(y, yb);

        arf_mul_ui(z2, z, 1025, MAG_BITS, ARF_RND_UP);
        arf_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(z, y) <= 0 && arf_cmpabs(y, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); arf_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_neg_log(xb, xb);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(z2);

        mag_clear(xb);
        mag_clear(yb);

        fmpz_clear(n);
    }

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, z2;
        fmpz_t n;
        mag_t xb, yb;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(z2);

        mag_init(xb);
        mag_init(yb);

        fmpz_init(n);

        mag_randtest_special(xb, state, 25);
        mag_randtest_special(yb, state, 25);

        if (n_randint(state, 2))
        {
            mag_neg_log_lower(yb, xb);

            arf_set_mag(x, xb);
            if (mag_cmp_2exp_si(xb, 0) >= 0)
                arf_zero(z);
            else
            {
                arf_log(z, x, MAG_BITS, ARF_RND_DOWN);
                arf_neg(z, z);
            }
        }
        else
        {
            arf_set_mag(x, xb);

            fmpz_randtest(n, state, 100);
            mag_mul_2exp_fmpz(xb, xb, n);
            mag_neg_log_lower(yb, xb);

            if (mag_cmp_2exp_si(xb, 0) >= 0)
            {
                arf_zero(z);
            }
            else
            {
                arf_log(z, x, 2 * MAG_BITS, ARF_RND_DOWN);
                arf_neg(z, z);
                arf_set_ui(z2, 2);
                arf_log(z2, z2, 2 * MAG_BITS, ARF_RND_DOWN);
                arf_submul_fmpz(z, z2, n, 2 * MAG_BITS, ARF_RND_DOWN);
            }
        }

        arf_set_mag(y, yb);

        arf_mul_ui(z2, z, 1023, MAG_BITS, ARF_RND_DOWN);
        arf_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(z2, y) <= 0 && arf_cmpabs(y, z) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); arf_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_neg_log_lower(xb, xb);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(z2);

        mag_clear(xb);
        mag_clear(yb);

        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
