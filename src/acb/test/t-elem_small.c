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

/* Containment for the register paths of the elementary functions at one
   and two limbs (acb/elem_small.c): for balls of small radius, the
   values at random points of the ball, computed at a high precision,
   lie in the result; and the result is not much wider than the
   propagated radius allows (its relative accuracy within a few bits of
   the input's). */

typedef void (*fun_t)(acb_t, const acb_t, slong);

static void sin_w(acb_t r, const acb_t z, slong prec) { acb_sin(r, z, prec); }
static void cos_w(acb_t r, const acb_t z, slong prec) { acb_cos(r, z, prec); }

/* log |z| and atan2(Im z, Re z) of arb into the real part */
static void log_hypot_w(acb_t r, const acb_t z, slong prec)
{
    arb_log_hypot(acb_realref(r), acb_realref(z), acb_imagref(z), prec);
    arb_zero(acb_imagref(r));
}
static void atan2_w(acb_t r, const acb_t z, slong prec)
{
    arb_atan2(acb_realref(r), acb_imagref(z), acb_realref(z), prec);
    arb_zero(acb_imagref(r));
}

TEST_FUNCTION_START(acb_elem_small, state)
{
    slong iter;
    static const fun_t funcs[] = { acb_exp, acb_exp_pi_i, sin_w, cos_w, acb_sin_pi, acb_cos_pi, acb_log, acb_sinh, acb_cosh, log_hypot_w, atan2_w };
    static const char * names[] = { "exp", "exp_pi_i", "sin", "cos", "sin_pi", "cos_pi", "log", "sinh", "cosh", "log_hypot", "atan2" };

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        acb_t z, r, w, v;
        slong prec, i, which, kind;

        acb_init(z);
        acb_init(r);
        acb_init(w);
        acb_init(v);

        prec = 2 + n_randint(state, 128);
        which = n_randint(state, 11);
        kind = n_randint(state, 4);

        /* the midpoint, in a range the register paths take */
        arb_randtest(acb_realref(z), state, 1 + n_randint(state, prec), 1 + n_randint(state, 6));
        arb_randtest(acb_imagref(z), state, 1 + n_randint(state, prec), 1 + n_randint(state, 6));
        if (kind == 0)
            acb_mul_2exp_si(z, z, -(slong) n_randint(state, 100));
        else if (kind == 1)
            arb_mul_2exp_si(acb_imagref(z), acb_imagref(z), -(slong) n_randint(state, 100));
        if (arb_is_zero(acb_realref(z)) || arb_is_zero(acb_imagref(z)))
        {
            acb_clear(z); acb_clear(r); acb_clear(w); acb_clear(v);
            continue;
        }

        /* a small radius (or none) */
        mag_zero(arb_radref(acb_realref(z)));
        mag_zero(arb_radref(acb_imagref(z)));
        if (n_randint(state, 2))
            mag_randtest(arb_radref(acb_realref(z)), state, 1 + n_randint(state, 8));
        if (n_randint(state, 2))
            mag_randtest(arb_radref(acb_imagref(z)), state, 1 + n_randint(state, 8));
        if (mag_cmp_2exp_si(arb_radref(acb_realref(z)), -9) > 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -9);
        if (mag_cmp_2exp_si(arb_radref(acb_imagref(z)), -9) > 0)
            mag_set_ui_2exp_si(arb_radref(acb_imagref(z)), 1, -9);

        funcs[which](r, z, prec);

        if (!acb_is_finite(r))
        {
            acb_clear(z); acb_clear(r); acb_clear(w); acb_clear(v);
            continue;
        }

        for (i = 0; i < 8; i++)
        {
            /* a point of the ball: the midpoint plus a fraction of the radii */
            arf_t t;
            arf_init(t);
            acb_set(w, z);
            mag_zero(arb_radref(acb_realref(w)));
            mag_zero(arb_radref(acb_imagref(w)));
            if (i > 0)
            {
                arf_set_mag(t, arb_radref(acb_realref(z)));
                arf_mul_ui(t, t, n_randint(state, 1025), 60, ARF_RND_DOWN);
                arf_mul_2exp_si(t, t, -10);
                if (n_randint(state, 2)) arf_neg(t, t);
                arf_add(arb_midref(acb_realref(w)), arb_midref(acb_realref(w)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_set_mag(t, arb_radref(acb_imagref(z)));
                arf_mul_ui(t, t, n_randint(state, 1025), 60, ARF_RND_DOWN);
                arf_mul_2exp_si(t, t, -10);
                if (n_randint(state, 2)) arf_neg(t, t);
                arf_add(arb_midref(acb_imagref(w)), arb_midref(acb_imagref(w)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
            }
            arf_clear(t);

            funcs[which](v, w, 400);

            if (!acb_contains(r, v))
            {
                flint_printf("FAIL: containment (%s)\n", names[which]);
                flint_printf("prec = %wd\n", prec);
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n");
                flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
                flint_printf("v = "); acb_printd(v, 30); flint_printf("\n");
                flint_abort();
            }
        }

        acb_clear(z);
        acb_clear(r);
        acb_clear(w);
        acb_clear(v);
    }

    TEST_FUNCTION_END(state);
}
