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
#include "mp_real.h"
#include "mp_real/impl.h"

/* Every run-time entry of the dense tables (dense_rt.c), in both tiers,
   against arb: a lower bound for log(1 + x), -log(1 - x) resp. atan(x)
   within one ulp.  (The kernels absorb a few ulps, so a wrong entry
   shows in their results only as an occasional error of a unit in a
   middle limb -- as with the last-level log1p entry of argument
   2^-64 - 2^-128, whose leading limb is zero, once computed with an
   undefined clz of zero.) */

TEST_FUNCTION_START(mp_real_dense_rt_tables, state)
{
#if FLINT_BITS == 64
    int which, tier;
    slong k, j, w, i;
    ulong r[MP_REAL_DENSE_RT_MAX + 1], a[2];
    arb_t x, f, d;
    fmpz_t t;

    arb_init(x);
    arb_init(f);
    arb_init(d);
    fmpz_init(t);

    for (tier = 0; tier < 2; tier++)
    {
        /* the tiers serve n <= 21 resp. beyond */
        slong n = tier ? MP_REAL_DENSE_RT_MAX : MP_REAL_DENSE_W + 1;

        for (which = 0; which < 3; which++)
        {
            for (k = 3; k <= MP_REAL_DENSE_RT_LEVELS; k++)
            {
                for (j = 0; j <= _mp_real_dense_rt_jmax(which, k); j++)
                {
                    slong prec;

                    w = _mp_real_dense_rt_entry(r, a, which, n, k, j);
                    prec = FLINT_BITS * w + 128;

                    /* x = (a1 B + a0) / B^2 exactly */
                    fmpz_set_uiui(t, a[1], a[0]);
                    arb_set_fmpz(x, t);
                    arb_mul_2exp_si(x, x, -2 * FLINT_BITS);

                    if (which == MP_REAL_DENSE_RT_EXP)
                        arb_log1p(f, x, prec);
                    else if (which == MP_REAL_DENSE_RT_LOG1P)
                    {
                        arb_neg(f, x);
                        arb_log1p(f, f, prec);
                        arb_neg(f, f);
                    }
                    else
                        arb_atan(f, x, prec);

                    /* d = f B^w - entry, required in [0, 1] */
                    fmpz_set_ui_array(t, r, w);
                    arb_mul_2exp_si(f, f, FLINT_BITS * w);
                    arb_sub_fmpz(d, f, t, prec);

                    if (!(arf_sgn(arb_midref(d)) >= 0 && arf_cmp_si(arb_midref(d), 1) <= 0)
                            || mag_cmp_2exp_si(arb_radref(d), -60) > 0)
                    {
                        flint_printf("FAIL: family %d, tier %d, k = %wd, j = %wd, "
                            "x = (%wu, %wu): f B^w - entry = ", which, tier, k, j,
                            a[1], a[0]);
                        arb_printd(d, 10);
                        flint_printf("\n");
                        for (i = w - 1; i >= 0; i--)
                            flint_printf(" %016wx", r[i]);
                        flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }
    }

    arb_clear(x);
    arb_clear(f);
    arb_clear(d);
    fmpz_clear(t);
#endif

    TEST_FUNCTION_END(state);
}
