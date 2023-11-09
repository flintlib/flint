/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

static mp_limb_t
refimpl(slong * exp, const fmpz_t x, int bits)
{
    fmpz_t t;
    slong xbits;
    mp_limb_t m;

    xbits = fmpz_bits(x);

    fmpz_init(t);
    fmpz_abs(t, x);

    if (xbits >= bits)
        fmpz_tdiv_q_2exp(t, t, xbits - bits);
    else
        fmpz_mul_2exp(t, t, bits - xbits);

    m = fmpz_get_ui(t);
    fmpz_clear(t);

    *exp = xbits - bits;

    return m;
}

TEST_FUNCTION_START(fmpz_abs_lbound_ui_2exp, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_t x;
        slong bits;
        slong exp, yexp;
        mp_limb_t yman, man;

        fmpz_init(x);
        fmpz_randtest_not_zero(x, state, 1 + n_randint(state, 400));

        bits = 1 + n_randint(state, FLINT_BITS - 1);

        yman = refimpl(&yexp, x, bits);
        man = fmpz_abs_lbound_ui_2exp(&exp, x, bits);

        if (FLINT_BIT_COUNT(man) != bits || (man != yman) || (exp != yexp))
        {
            flint_printf("FAIL\n");
            flint_printf("bits = %wd, count = %u\n\n", bits, FLINT_BIT_COUNT(man));
            flint_printf("x = "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("bits(x) = %wd\n\n", fmpz_bits(x));
            flint_printf("man = %wu, exp = %wd\n", man, exp);
            flint_printf("yman = %wu, yexp = %wd\n", yman, yexp);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}
