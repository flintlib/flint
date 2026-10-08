/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"

/* mp_real_exp_notab_log2, mp_real_exp_notab_squaring and mp_real_exp_agm
   against arb_exp: the outputs contain exp of the input ball, and for
   exact inputs have a relative accuracy of about n limbs, for arguments
   from tiny to 2^31 (2^26 on 32-bit machines) in absolute value of
   either sign, at precisions
   where the AGM step is taken (from about 30 limbs) and where it falls
   back to the squarings. */

static void
_exp_huge_set_arf(mp_real_t x, const arf_t m)
{
    fmpz_t man, e;

    if (arf_is_zero(m))
    {
        mp_real_zero(x);
        return;
    }
    fmpz_init(man);
    fmpz_init(e);
    arf_get_fmpz_2exp(man, e, m);
    mp_real_set_fmpz(x, man);
    mp_real_mul_2exp_si(x, x, fmpz_get_si(e));
    fmpz_clear(man);
    fmpz_clear(e);
}

TEST_FUNCTION_START(mp_real_exp_huge, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong n, i;
        int kind, method;
        arf_t m;
        mp_real_t x, y;
        arb_t a, r, ref;

        n = 1 + n_randint(state, (iter % 30 == 0) ? 600 : (iter % 3 == 0) ? 100 : 20);
        kind = n_randint(state, 5);

        arf_init(m);
        mp_real_init(x);
        mp_real_init(y);
        arb_init(a);
        arb_init(r);
        arb_init(ref);

        switch (kind)
        {
            case 0:
                arf_randtest(m, state, FLINT_BITS * n + 10, 4);
                break;
            case 1:     /* |m| < 2^31, within the safe range 2^(FLINT_BITS - 5) */
                arf_randtest(m, state, FLINT_BITS * n + 10, 5);
                if (!arf_is_zero(m) && arf_cmpabs_2exp_si(m, FLINT_BITS - 6) >= 0)
                    arf_mul_2exp_si(m, m, FLINT_BITS - 6 - fmpz_get_si(ARF_EXPREF(m)));
                break;
            case 2:     /* tiny */
                arf_randtest(m, state, FLINT_BITS * n + 10, 3);
                arf_mul_2exp_si(m, m, -(slong) n_randint(state, FLINT_BITS * n + 100));
                break;
            case 3:
                arf_set_si(m, (slong) n_randint(state, 2000) - 1000);
                break;
            default:    /* large, short: as for the exp(C) of HRR */
                arf_set_ui(m, 1 + n_randint(state, 1000000));
                arf_mul_2exp_si(m, m, -(slong) n_randint(state, 4));
                if (n_randint(state, 2))
                    arf_neg(m, m);
                break;
        }

        _exp_huge_set_arf(x, m);
        if (x->size != 0 && n_randint(state, 4) == 0)
            x->err = 1 + n_randint(state, 100);
        mp_real_get_arb(a, x);
        arb_exp(ref, a, FLINT_BITS * n + 200);

        for (i = 0; i < 2; i++)
        {
            method = (i == 0) ? n_randint(state, 3) : n_randint(state, 3);
            if (method == 0)
                mp_real_exp_notab_log2(y, x, n);
            else if (method == 1)
                mp_real_exp_notab_squaring(y, x, n);
            else
                mp_real_exp_agm(y, x, n);
            mp_real_get_arb(r, y);

            if (!arb_overlaps(r, ref)
                || (x->err == 0 && arb_rel_accuracy_bits(r) < FLINT_BITS * n - 8))
            {
                flint_printf("FAIL: method %d, n = %wd, kind %d\n", method, n, kind);
                flint_printf("x = "); arb_printd(a, 30); flint_printf("\n");
                flint_printf("y = "); arb_printd(r, 30); flint_printf("\n");
                flint_printf("    "); arb_printd(ref, 30); flint_printf("\n");
                flint_abort();
            }
        }

        arf_clear(m);
        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(a);
        arb_clear(r);
        arb_clear(ref);
    }

    TEST_FUNCTION_END(state);
}
