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

/* mp_real_sin_cos_bits: the output balls contain sin and cos of a
   random point of the input ball, for tiny, unit, one-limb, many-limb,
   near k pi/2, and inexact arguments of either sign; for exact
   arguments the sine has a relative accuracy of about prec bits when
   |x| < 1 and both outputs an absolute accuracy of about prec bits;
   one output may be NULL, and an output may alias the input. */

static double
sc_abs_acc(const arb_t a)
{
    if (mag_is_zero(arb_radref(a)))
        return 1e9;
    return -mag_get_d_log2_approx(arb_radref(a));
}

TEST_FUNCTION_START(mp_real_sin_cos_bits, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        slong prec, L, e2, ip;
        int kind, inexact, mode;
        mp_real_t x, s, c, y;
        arb_t xa, sa, ca, t, ts, tc;
        arf_t pt;
        nn_ptr p;

        if (iter % 100 == 0)
            prec = 2 + n_randint(state, 12000);
        else if (iter % 10 == 0)
            prec = 2 + n_randint(state, 2500);
        else
            prec = 2 + n_randint(state, 300);

        kind = n_randint(state, 7);

        mp_real_init(x);
        mp_real_init(s);
        mp_real_init(c);
        mp_real_init(y);
        arb_init(xa);
        arb_init(sa);
        arb_init(ca);
        arb_init(t);
        arb_init(ts);
        arb_init(tc);
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
            case 2: e2 = 1 + n_randint(state, 64); break;
            case 3: e2 = 1 + n_randint(state, 3000); break;
            default: e2 = (slong) n_randint(state, 200) - 100; break;
        }
        _mp_real_set_mpn_2exp(x, p, L, e2 - FLINT_BITS * L);

        if (kind == 4)
        {
            /* near k pi/2 */
            mp_real_t d;
            mp_real_init(d);
            mp_real_const_pi4(x, L + 4, 1);
            mp_real_mul_2exp_si(x, x, 1);
            mp_real_mul_ui(x, x, n_randint(state, 1000000) + 1, L + 4);
            x->err = 0;
            if (n_randint(state, 2))
            {
                mp_real_set_ui(d, 1);
                mp_real_mul_2exp_si(d, d, -(slong) n_randint(state, 200));
                mp_real_add(x, x, d, L + 10);
                x->err = 0;
            }
            mp_real_clear(d);
        }

        if (n_randint(state, 2))
            mp_real_neg(x, x);

        inexact = (n_randint(state, 4) == 0);
        if (inexact)
            mp_real_add_rel_error_2exp_si(x, -(slong) n_randint(state, prec + 20));

        mp_real_get_arb(xa, x);

        /* 0: both outputs; 1, 2: one NULL; 3, 4: an output aliases x */
        mode = n_randint(state, 5);
        if (mode == 0)
            mp_real_sin_cos_bits(s, c, x, prec);
        else if (mode == 1)
        {
            mp_real_sin_cos_bits(s, NULL, x, prec);
            mp_real_sin_cos_bits(NULL, c, x, prec);
        }
        else if (mode == 2)
        {
            mp_real_sin_cos_bits(NULL, c, x, prec);
            mp_real_sin_cos_bits(s, NULL, x, prec);
        }
        else if (mode == 3)
        {
            mp_real_set(y, x);
            mp_real_sin_cos_bits(y, c, y, prec);
            mp_real_swap(s, y);
        }
        else
        {
            mp_real_set(y, x);
            mp_real_sin_cos_bits(s, y, y, prec);
            mp_real_swap(c, y);
        }

        mp_real_get_arb(sa, s);
        mp_real_get_arb(ca, c);

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
        ip = arf_is_zero(pt) ? 0 : FLINT_MAX(0, fmpz_get_si(ARF_EXPREF(pt)));
        arb_sin_cos(ts, tc, t, prec + 100 + ip);

        if (!arb_overlaps(ts, sa) || !arb_overlaps(tc, ca))
        {
            flint_printf("FAIL: containment, kind %d, prec %wd, inexact %d, mode %d\n",
                kind, prec, inexact, mode);
            flint_printf("x = "); arb_printd(xa, 30);
            flint_printf("\nsin = "); arb_printd(sa, 30);
            flint_printf("\n      "); arb_printd(ts, 30);
            flint_printf("\ncos = "); arb_printd(ca, 30);
            flint_printf("\n      "); arb_printd(tc, 30);
            flint_printf("\n");
            flint_abort();
        }

        if (!inexact)
        {
            if (arf_cmpabs_2exp_si(arb_midref(xa), 0) < 0
                && arb_rel_accuracy_bits(sa) < prec - 4)
            {
                flint_printf("FAIL: sin relative accuracy, kind %d, prec %wd, acc %wd\n",
                    kind, prec, arb_rel_accuracy_bits(sa));
                flint_printf("x = "); arb_printd(xa, 30); flint_printf("\n");
                flint_abort();
            }

            if (FLINT_MIN(sc_abs_acc(sa), sc_abs_acc(ca)) < prec - 4)
            {
                flint_printf("FAIL: absolute accuracy, kind %d, prec %wd\n", kind, prec);
                flint_printf("x = "); arb_printd(xa, 30);
                flint_printf("\nsin = "); arb_printd(sa, 30);
                flint_printf("\ncos = "); arb_printd(ca, 30);
                flint_printf("\n");
                flint_abort();
            }
        }

        flint_free(p);
        mp_real_clear(x);
        mp_real_clear(s);
        mp_real_clear(c);
        mp_real_clear(y);
        arb_clear(xa);
        arb_clear(sa);
        arb_clear(ca);
        arb_clear(t);
        arb_clear(ts);
        arb_clear(tc);
        arf_clear(pt);
    }

    /* huge arguments and zero */
    {
        mp_real_t x, s, c;
        arb_t a;

        mp_real_init(x);
        mp_real_init(s);
        mp_real_init(c);
        arb_init(a);

        mp_real_set_ui(x, 3);
        mp_real_mul_2exp_si(x, x, 70000);
        mp_real_sin_cos_bits(s, c, x, 64);
        mp_real_get_arb(a, s);
        if (!arb_contains_si(a, 1) || !arb_contains_si(a, -1))
            TEST_FUNCTION_FAIL("huge argument\n");

        mp_real_zero(x);
        mp_real_sin_cos_bits(s, c, x, 64);
        mp_real_get_arb(a, c);
        if (!mp_real_is_zero(s) || !arb_is_one(a))
            TEST_FUNCTION_FAIL("zero\n");

        mp_real_clear(x);
        mp_real_clear(s);
        mp_real_clear(c);
        arb_clear(a);
    }

    TEST_FUNCTION_END(state);
}
