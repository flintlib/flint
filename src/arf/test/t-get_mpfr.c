/*
    Copyright (C) 2012, 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arf.h"

TEST_FUNCTION_START(arf_get_mpfr, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits;
        arf_t x, z;
        mpfr_t y;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(z);
        mpfr_init2(y, bits);

        arf_randtest_special(x, state, bits, 10);
        arf_get_mpfr(y, x, MPFR_RNDN);
        arf_set_mpfr(z, y);

        if (!arf_equal(x, z))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(z);
        mpfr_clear(y);
    }

    /* test set_mpfr out of range */
    {
        arf_t x, z;
        mpfr_t y;
        fmpz_t e;
        int r;

        fmpz_init(e);
        arf_init(x);
        arf_init(z);
        mpfr_init2(y, 53);
        fmpz_one(e);
        fmpz_mul_2exp(e, e, 100);

        arf_set_si(x, 1);
        arf_mul_2exp_fmpz(x, x, e);
        r = arf_get_mpfr(y, x, MPFR_RNDN);
        arf_set_mpfr(x, y);

        if (!mpfr_inf_p(y) || mpfr_sgn(y) <= 0 || !mpfr_overflow_p())
        {
            flint_printf("expected +inf with overflow\n\n");
            arf_print(z); flint_printf("\n\n");
            flint_printf("r = %d \n\n", r);
            flint_abort();
        }

        mpfr_clear_flags();

        arf_set_si(x, -1);
        arf_mul_2exp_fmpz(x, x, e);
        r = arf_get_mpfr(y, x, MPFR_RNDN);
        arf_set_mpfr(x, y);

        if (!mpfr_inf_p(y) || mpfr_sgn(y) >= 0 || !mpfr_overflow_p())
        {
            flint_printf("expected -inf with overflow\n\n");
            arf_print(z); flint_printf("\n\n");
            flint_printf("r = %d \n\n", r);
            flint_abort();
        }

        mpfr_clear_flags();

        fmpz_neg(e, e);

        arf_set_si(x, 1);
        arf_mul_2exp_fmpz(x, x, e);
        r = arf_get_mpfr(y, x, MPFR_RNDN);
        arf_set_mpfr(x, y);

        if (!mpfr_zero_p(y) || mpfr_signbit(y) || !mpfr_underflow_p())
        {
            flint_printf("expected +0 with underflow\n\n");
            arf_print(z); flint_printf("\n\n");
            flint_printf("r = %d \n\n", r);
            flint_abort();
        }

        mpfr_clear_flags();

        arf_set_si(x, -1);
        arf_mul_2exp_fmpz(x, x, e);
        r = arf_get_mpfr(y, x, MPFR_RNDN);
        arf_set_mpfr(x, y);

        if (!mpfr_zero_p(y) || !mpfr_signbit(y) || !mpfr_underflow_p())
        {
            flint_printf("expected -0 with underflow\n\n");
            arf_print(z); flint_printf("\n\n");
            flint_printf("r = %d \n\n", r);
            flint_abort();
        }

        mpfr_clear_flags();

        arf_clear(x);
        arf_clear(z);
        mpfr_clear(y);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
