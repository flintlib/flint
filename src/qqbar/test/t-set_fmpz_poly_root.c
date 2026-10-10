/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "fmpz_poly.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_set_fmpz_poly_root, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        acb_t z;
        slong prec;
        int ok;

        qqbar_init(x);
        qqbar_init(y);
        acb_init(z);

        qqbar_randtest(x, state, 6, 10);
        prec = 10 + n_randint(state, 200);
        qqbar_get_acb(z, x, prec);

        /* a wide enclosure around a conjugate may or may not isolate it */
        if (n_randint(state, 2))
        {
            arb_add_error_2exp_si(acb_realref(z), -n_randint(state, prec));
            arb_add_error_2exp_si(acb_imagref(z), -n_randint(state, prec));
        }

        ok = qqbar_set_fmpz_poly_root(y, QQBAR_POLY(x), z, 2 * prec);

        if (ok && !qqbar_equal(x, y))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 20); flint_printf("\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n");
            flint_abort();
        }

        /* a square enclosure much narrower than the root separation
           (beyond the Mahler bound 2^-(d (h + 4))) with the root well
           inside is certified. (The rounded enclosures of qqbar_get_acb
           need not qualify: the root may lie on the boundary, as for
           x^2 - 29x + 629 at 14.4998779296875 +/- 2^-13, the enclosure
           may be exact (i), or much narrower in one direction than in
           the other, which interval Newton on rectangles does not handle.) */
        if (qqbar_degree(x) >= 2)
        {
            slong hprec = 2 * qqbar_degree(x) * (qqbar_height_bits(x) + 4) + 64;

            qqbar_get_acb(z, x, 2 * hprec);
            acb_get_mid(z, z);
            if (!qqbar_is_real(x))
                arb_add_error_2exp_si(acb_imagref(z), -hprec);
            if (!arb_is_zero(acb_realref(z)) || qqbar_is_real(x))
                arb_add_error_2exp_si(acb_realref(z), -hprec);

            ok = qqbar_set_fmpz_poly_root(y, QQBAR_POLY(x), z, 2 * hprec);

            if (!ok || !qqbar_equal(x, y))
            {
                flint_printf("FAIL (accurate enclosure not certified)\n\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n");
                flint_printf("z = "); acb_printd(z, 20); flint_printf("\n");
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(y);
        acb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
