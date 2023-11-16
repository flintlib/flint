/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "mag.h"
#include "arb.h"

TEST_FUNCTION_START(mag_hurwitz_zeta_uiui, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mag_t zb;
        slong prec, k;
        ulong s, a;
        arb_t t, r, z;

        mag_init(zb);
        arb_init(t);
        arb_init(r);
        arb_init(z);

        mag_randtest_special(zb, state, 6);
        s = n_randtest(state);
        a = n_randtest(state);
        prec = MAG_BITS + 2 * FLINT_BIT_COUNT(s);

        mag_hurwitz_zeta_uiui(zb, s, a);

        arb_zero(r);
        for (k = 0; k < 50; k++)
        {
            arb_set_ui(t, a);
            arb_add_ui(t, t, k, prec);
            arb_pow_ui(t, t, s, prec);
            arb_inv(t, t, prec);
            arb_add(r, r, t, prec);
        }

        arb_zero(z);
        mag_set(arb_radref(z), zb);
        if (!arb_is_finite(z))
            arb_indeterminate(z);

        if (!arb_contains(z, r))
        {
            flint_printf("FAIL\n\n");
            flint_printf("s = %wu\n\n", s);
            flint_printf("a = %wd\n\n", a);
            flint_printf("zb = "); mag_printd(zb, 15); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 15); flint_printf("\n\n");
            flint_printf("r = "); arb_printd(r, 15); flint_printf("\n\n");
            flint_abort();
        }

        mag_clear(zb);
        arb_clear(t);
        arb_clear(r);
        arb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
