/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "arb/elem.h"

/* _arb_set_mpn_err and the register frames _arb_small_set_1/2: the ball
   contains the frame Y 2^e plus or minus the error count err (as an
   exact ball), for frames with leading zero limbs and bits (down to a
   zero frame) and trailing zero limbs and bits, errors from 0 to
   beyond 2^20, and precisions on both sides of the frame's size; the
   register versions have a radius within a factor 2 of that of
   _arb_set_mpn_err plus an ulp. */
TEST_FUNCTION_START(arb_set_mpn_err, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong y[5], t[5];
        slong len, i, e, prec;
        ulong err;
        int sgnbit;
        arb_t res, ref, exact;
        fmpz_t Y;
        mag_t m;

        arb_init(res);
        arb_init(ref);
        arb_init(exact);
        fmpz_init(Y);
        mag_init(m);

        len = 3 + n_randint(state, 2);
        for (i = 0; i < len; i++)
            y[i] = n_randtest(state);
        /* leading zero limbs and bits */
        for (i = len - 1; i >= 0 && n_randint(state, 3) == 0; i--)
            y[i] = 0;
        if (i >= 0 && n_randint(state, 2))
            y[i] >>= n_randint(state, FLINT_BITS);
        /* trailing zero limbs and bits (exact and nearly exact frames) */
        if (n_randint(state, 3) == 0)
        {
            slong j;
            for (j = 0; j < len - 1 && n_randint(state, 3) != 0; j++)
                y[j] = 0;
            if (n_randint(state, 2))
                y[j] &= ~((UWORD(1) << n_randint(state, FLINT_BITS)) - 1);
        }

        e = (slong) n_randint(state, 400) - 200;
        err = n_randint(state, 4) == 0 ? 0 : n_randtest(state) >> n_randint(state, FLINT_BITS);
        if (n_randint(state, 2))
            err &= (UWORD(1) << 19) - 1;
        sgnbit = n_randint(state, 2);
        prec = 2 + n_randint(state, 2 * FLINT_BITS + 10);

        /* the exact frame plus or minus err, in units of 2^e */
        fmpz_set_ui_array(Y, y, len);
        arb_set_fmpz(exact, Y);
        if (sgnbit)
            arb_neg(exact, exact);
        mag_set_ui(arb_radref(exact), err);
        arb_mul_2exp_si(exact, exact, e);

        for (i = 0; i < len; i++)
            t[i] = y[i];
        _arb_set_mpn_err(ref, t, len, e, err, err, sgnbit, prec);
        if (!arb_contains(ref, exact))
        {
            flint_printf("FAIL: _arb_set_mpn_err (containment)\n");
            arb_printd(ref, 30); flint_printf("\n");
            arb_printd(exact, 30); flint_printf("\n");
            flint_abort();
        }

#if FLINT_BITS == 64
        if (len == 3 || len == 4)
        {
            if (len == 3)
                _arb_small_set_1(res, y[2], y[1], y[0], e + 2 * FLINT_BITS, err, sgnbit, prec);
            else
                _arb_small_set_2(res, y[3], y[2], y[1], y[0], e + 3 * FLINT_BITS, err, sgnbit, prec);

            mag_mul_2exp_si(m, arb_radref(ref), 1);
            if (!arf_is_special(arb_midref(ref)))
                arf_mag_add_ulp(m, m, arb_midref(ref), prec);

            if (!arb_contains(res, exact) || mag_cmp(arb_radref(res), m) > 0)
            {
                flint_printf("FAIL: _arb_small_set_%wd, prec = %wd, err = %wu\n", len - 2, prec, err);
                arb_printd(res, 30); flint_printf("\n");
                arb_printd(ref, 30); flint_printf("\n");
                arb_printd(exact, 30); flint_printf("\n");
                flint_abort();
            }
        }
#endif

        arb_clear(res);
        arb_clear(ref);
        arb_clear(exact);
        fmpz_clear(Y);
        mag_clear(m);
    }

#if FLINT_BITS == 64
    /* the frame (2^63, 1, 0), a zero low limb: 2^63 + 2^-64 for k = 0,
       inexact below 128 bits (with prec = 64 the rounding point at the
       limb boundary) and exact from 128 bits, with a zero radius for
       err = 0 */
    {
        arb_t res, exact;
        slong prec;
        ulong err;
        arb_init(res);
        arb_init(exact);
        for (prec = 40; prec <= 128; prec += 8)
        {
            for (err = 0; err <= 3; err += 3)
            {
                /* 2^63 + 2^-64 plus or minus err 2^-128 */
                arb_set_ui(exact, 1);
                arb_mul_2exp_si(exact, exact, -127);
                arb_add_ui(exact, exact, 1, ARF_PREC_EXACT);
                arb_mul_2exp_si(exact, exact, 63);
                mag_set_ui_2exp_si(arb_radref(exact), err, -128);

                _arb_small_set_1(res, UWORD(1) << 63, UWORD(1), 0, 0, err, 0, prec);
                if (!arb_contains(res, exact) || (prec == 128 && !arb_equal(res, exact)))
                {
                    flint_printf("FAIL: frame with a zero low limb, prec = %wd, err = %wu\n", prec, err);
                    arb_printd(res, 40); flint_printf("\n");
                    arb_printd(exact, 40); flint_printf("\n");
                    flint_abort();
                }
            }
        }
        arb_clear(res);
        arb_clear(exact);
    }
#endif

    TEST_FUNCTION_END(state);
}
