/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"

/* mp_real_exp_bits, mp_real_log_bits, mp_real_atan_bits: the output
   ball contains f at a random point of the input ball, for tiny, unit,
   one-limb and many-limb arguments (up to 2^49 for exp), arguments
   near 1 (and near -1 for atan), small integers and powers of two,
   inexact arguments and aliased outputs; exact arguments give a
   relative accuracy of about prec bits; log reports non-positive balls */

TEST_FUNCTION_START(mp_real_exp_log_atan_bits, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        int func = n_randint(state, 3), kind, inexact, alias, ok = 1;
        slong prec, L, e2, wp;
        mp_real_t x, y, d;
        arb_t xa, ya, t, ty;
        arf_t pt;
        nn_ptr p;

        if (iter % 100 == 0)
            prec = 2 + n_randint(state, 12000);
        else if (iter % 10 == 0)
            prec = 2 + n_randint(state, 2500);
        else
            prec = 2 + n_randint(state, 300);

        kind = n_randint(state, 8);

        mp_real_init(x);
        mp_real_init(y);
        mp_real_init(d);
        arb_init(xa);
        arb_init(ya);
        arb_init(t);
        arb_init(ty);
        arf_init(pt);

        L = 1 + n_randint(state, prec / FLINT_BITS + 3);
        p = flint_malloc(L * sizeof(ulong));
        flint_mpn_rrandom(p, state, L);
        if (p[L - 1] == 0)
            p[L - 1] = 1;

        switch (kind)
        {
            case 0: e2 = -(slong) n_randint(state, prec + 200); break;
            case 1: e2 = -(slong) n_randint(state, 64); break;
            case 2: e2 = 1 + n_randint(state, func == 0 ? 20 : 64); break;
            case 3: e2 = (func == 0) ? (slong) n_randint(state, FLINT_BITS - 14)
                        : (slong) n_randint(state, 3000) - 1500; break;
            default: e2 = (slong) n_randint(state, 20) - 10; break;
        }
        _mp_real_set_mpn_2exp(x, p, L, e2 - FLINT_BITS * L);

        if (kind == 4 || kind == 5)
        {
            /* 1 +- 2^-k (+ a perturbation below) */
            mp_real_set_ui(d, 1);
            mp_real_mul_2exp_si(d, d, -(slong) n_randint(state,
                (kind == 4) ? prec + 100 : 64) - 1);
            mp_real_set(y, x);
            mp_real_set_ui(x, 1);
            if (n_randint(state, 2))
                mp_real_add(x, x, d, 10000);
            else
                mp_real_sub(x, x, d, 10000);
            if (n_randint(state, 2))
            {
                mp_real_mul_2exp_si(y, y, -(slong) n_randint(state, prec + 300) - 2);
                mp_real_add(x, x, y, 10000);
            }
            x->err = 0;
        }

        if (kind == 6 && func != 1)
        {
            mp_real_set_ui(x, n_randint(state, 5));
            mp_real_mul_2exp_si(x, x, (slong) n_randint(state, 20) - 10);
        }

        if (func != 1 && n_randint(state, 2))
            mp_real_neg(x, x);
        if (func == 1 && n_randint(state, 50) == 0)
            mp_real_neg(x, x);

        inexact = (n_randint(state, 4) == 0);
        if (inexact)
        {
            if (mp_real_is_zero(x))
                mp_real_add_error_2exp_si(x, -(slong) n_randint(state, 100));
            else
                mp_real_add_error_2exp_si(x, mp_real_abs_bound_lt_2exp_si(x)
                    - 1 - (slong) n_randint(state, prec + 20));
        }
        mp_real_get_arb(xa, x);

        alias = (n_randint(state, 4) == 0);
        if (alias)
        {
            mp_real_set(y, x);
            if (func == 0)
                mp_real_exp_bits(y, y, prec);
            else if (func == 1)
                ok = mp_real_log_bits(y, y, prec);
            else
                mp_real_atan_bits(y, y, prec);
        }
        else
        {
            if (func == 0)
                mp_real_exp_bits(y, x, prec);
            else if (func == 1)
                ok = mp_real_log_bits(y, x, prec);
            else
                mp_real_atan_bits(y, x, prec);
        }
        mp_real_get_arb(ya, y);

        if (func == 1 && ok != arb_is_positive(xa))
        {
            flint_printf("FAIL: log positivity, ok = %d\nx = ", ok);
            arb_printd(xa, 30);
            flint_printf("\n");
            flint_abort();
        }

        if (ok)
        {
            /* a random point of the input ball */
            arf_set(pt, arb_midref(xa));
            if (inexact)
            {
                arf_t r;
                arf_init(r);
                arf_set_mag(r, arb_radref(xa));
                arf_mul_2exp_si(r, r, -1);
                if (n_randint(state, 2))
                    arf_neg(r, r);
                arf_add(pt, pt, r, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_clear(r);
            }
            arb_set_arf(t, pt);
            wp = prec + 100 + (arf_is_zero(pt) ? 0
                : FLINT_ABS(fmpz_get_si(ARF_EXPREF(pt))));
            if (func == 0)
                arb_exp(ty, t, wp);
            else if (func == 1)
                arb_log(ty, t, wp);
            else
                arb_atan(ty, t, wp);

            if (!arb_overlaps(ty, ya))
            {
                flint_printf("FAIL: containment, func %d, kind %d, prec %wd, "
                    "inexact %d, alias %d\nx = ", func, kind, prec, inexact, alias);
                arb_printd(xa, 30);
                flint_printf("\ny = "); arb_printd(ya, 30);
                flint_printf("\n    "); arb_printd(ty, 30);
                flint_printf("\n");
                flint_abort();
            }

            if (!inexact && !arb_is_exact(ya)
                && arb_rel_accuracy_bits(ya) < prec - 4)
            {
                flint_printf("FAIL: accuracy, func %d, kind %d, prec %wd, acc %wd\nx = ",
                    func, kind, prec, arb_rel_accuracy_bits(ya));
                arb_printd(xa, 30);
                flint_printf("\ny = "); arb_printd(ya, 30);
                flint_printf("\n");
                flint_abort();
            }
        }

        flint_free(p);
        mp_real_clear(x);
        mp_real_clear(y);
        mp_real_clear(d);
        arb_clear(xa);
        arb_clear(ya);
        arb_clear(t);
        arb_clear(ty);
        arf_clear(pt);
    }

    TEST_FUNCTION_END(state);
}
