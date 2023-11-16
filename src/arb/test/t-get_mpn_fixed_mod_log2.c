/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_get_mpn_fixed_mod_log2, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x;
        fmpz_t q;
        mp_ptr w;
        arb_t wb, t;
        mp_size_t wn;
        slong prec, prec2;
        int success;
        mp_limb_t error;

        prec = 2 + n_randint(state, 10000);
        wn = 1 + n_randint(state, 200);
        prec2 = FLINT_MAX(prec, wn * FLINT_BITS) + 100;

        arf_init(x);
        arb_init(wb);
        arb_init(t);
        fmpz_init(q);
        w = flint_malloc(sizeof(mp_limb_t) * wn);

        arf_randtest(x, state, prec, 14);

        /* this should generate numbers close to multiples of log(2) */
        if (n_randint(state, 4) == 0)
        {
            arb_const_log2(t, prec);
            fmpz_randtest(q, state, 200);
            arb_mul_fmpz(t, t, q, prec);
            arf_add(x, x, arb_midref(t), prec, ARF_RND_DOWN);
        }

        success = _arb_get_mpn_fixed_mod_log2(w, q, &error, x, wn);

        if (success)
        {
            _arf_set_mpn_fixed(arb_midref(wb), w, wn, wn, 0, FLINT_BITS * wn, ARB_RND);
            mag_set_ui_2exp_si(arb_radref(wb), error, -FLINT_BITS * wn);

            arb_const_log2(t, prec2);

            arb_mul_fmpz(t, t, q, prec2);
            arb_add(t, t, wb, prec2);

            if (!arb_contains_arf(t, x))
            {
                flint_printf("FAIL (containment)\n");
                flint_printf("x = "); arf_printd(x, 50); flint_printf("\n\n");
                flint_printf("q = "); fmpz_print(q); flint_printf("\n\n");
                flint_printf("w = "); arb_printd(wb, 50); flint_printf("\n\n");
                flint_printf("t = "); arb_printd(t, 50); flint_printf("\n\n");
                flint_abort();
            }

            arb_const_log2(t, prec2);

            if (arf_sgn(arb_midref(wb)) < 0 ||
                arf_cmp(arb_midref(wb), arb_midref(t)) >= 0)
            {
                flint_printf("FAIL (expected 0 <= w < log(2))\n");
                flint_printf("x = "); arf_printd(x, 50); flint_printf("\n\n");
                flint_printf("w = "); arb_printd(wb, 50); flint_printf("\n\n");
                flint_abort();
            }
        }

        flint_free(w);
        fmpz_clear(q);
        arf_clear(x);
        arb_clear(wb);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
