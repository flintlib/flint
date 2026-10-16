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
#include "ball_helpers.h"

/* mp_real_sin_cos_pi_bits and mp_real_tan_pi_bits: the outputs contain
   the values at a random point of the input ball (for tan when the
   status is 1); for exact arguments they have a relative accuracy of
   about prec bits at any magnitude and near every zero and pole (the
   arguments include k/4 exactly, k/2 + small, and huge and tiny values),
   the values 0, +-1 (and tan(pi/4) = 1) are exact, and the tangent's
   status is 0 only at a pole or for an inexact argument whose ball,
   widened three times, contains one; either sine/cosine output may be
   NULL, and an output may alias the input. */

/* x = m exactly */
static void
trig_pi_test_set_arf(mp_real_t x, const arf_t m)
{
    fmpz_t man, e;
    slong sz;
    nn_ptr d;
    int neg;

    if (arf_is_zero(m))
    {
        mp_real_zero(x);
        return;
    }

    fmpz_init(man);
    fmpz_init(e);
    arf_get_fmpz_2exp(man, e, m);
    neg = (fmpz_sgn(man) < 0);
    fmpz_abs(man, man);
    sz = fmpz_size(man);
    d = flint_malloc(sz * sizeof(ulong));
    fmpz_get_ui_array(d, sz, man);
    _mp_real_set_mpn_2exp(x, d, sz, fmpz_get_si(e));
    if (neg)
        mp_real_neg(x, x);
    flint_free(d);
    fmpz_clear(man);
    fmpz_clear(e);
}

/* the output got against the reference ref: containment, and for an
   exact argument exactness or the relative accuracy */
static void
trig_pi_test_check(const arb_t got, const arb_t ref, const arb_t xa,
    int inexact, slong prec, const char * what)
{
    if (!arb_overlaps(got, ref))
    {
        flint_printf("FAIL: %s containment, prec %wd\n", what, prec);
        flint_printf("x = "); arb_printd(xa, 30);
        flint_printf("\ny = "); arb_printd(got, 30);
        flint_printf("\n    "); arb_printd(ref, 30);
        flint_printf("\n");
        flint_abort();
    }

    if (!inexact && (arb_is_exact(ref) ? !arb_equal(got, ref)
            : arb_rel_accuracy_bits(got) < prec - 4))
    {
        flint_printf("FAIL: %s accuracy, prec %wd, acc %wd\n", what, prec,
            arb_rel_accuracy_bits(got));
        flint_printf("x = "); arb_printd(xa, 30);
        flint_printf("\ny = "); arb_printd(got, 30);
        flint_printf("\n    "); arb_printd(ref, 30);
        flint_printf("\n");
        flint_abort();
    }
}

TEST_FUNCTION_START(mp_real_trig_pi, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong prec, wp;
        int kind, inexact, alias, which;
        mp_real_t x, y1, y2;
        arb_t xa, ya, rs, rc, rt, t;
        arf_t m, pt;

        if (iter % 100 == 0)
            prec = 2 + n_randint(state, 12000);
        else if (iter % 10 == 0)
            prec = 2 + n_randint(state, 2500);
        else
            prec = 2 + n_randint(state, 300);

        kind = n_randint(state, 5);
        inexact = (n_randint(state, 3) == 0);
        alias = n_randint(state, 2);
        which = n_randint(state, 4);    /* both, sin, cos, tan */

        mp_real_init(x);
        mp_real_init(y1);
        mp_real_init(y2);
        arb_init(xa);
        arb_init(ya);
        arb_init(rs);
        arb_init(rc);
        arb_init(rt);
        arb_init(t);
        arf_init(m);
        arf_init(pt);

        switch (kind)
        {
            case 0:
                arf_randtest(m, state, prec + 10, 6);
                break;
            case 1:     /* huge or tiny */
                arf_randtest(m, state, prec + 10, 16);
                break;
            case 2:     /* k/4 exactly */
                arf_set_si(m, (slong) n_randint(state, 1000) - 500);
                arf_mul_2exp_si(m, m, -(slong) n_randint(state, 3));
                break;
            default:    /* k/2 + small */
            {
                slong d = n_randint(state, 3 * prec + 100);
                arf_t u;
                arf_init(u);
                arf_set_si(m, (slong) n_randint(state, 1000) - 500);
                arf_mul_2exp_si(m, m, -(slong) n_randint(state, 2));
                arf_set_ui_2exp_si(u, 1 + n_randint(state, 1000), -d);
                if (n_randint(state, 2))
                    arf_neg(u, u);
                arf_add(m, m, u, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_clear(u);
            }
        }

        trig_pi_test_set_arf(x, m);
        if (inexact)
        {
            slong e = (arf_is_zero(m) ? 0 : arf_abs_bound_lt_2exp_si(m))
                - (slong) n_randint(state, prec + 40);
            if (n_randint(state, 4) == 0)
                e = -(slong) n_randint(state, 5 * prec + 10);
            mp_real_add_error_2exp_si(x, e);
        }
        mp_real_get_arb(xa, x);

        mp_real_random_point(pt, x, state);
        wp = 2 * prec + 3 * FLINT_BITS * (x->size + 8) + 200;
        arb_set_arf(t, pt);
        arb_sin_cos_pi(rs, rc, t, wp);
        arb_tan_pi(rt, t, wp);

        if (which < 3)
        {
            mp_real_ptr os = (which != 2) ? y1 : NULL;
            mp_real_ptr oc = (which != 1) ? y2 : NULL;

            if (alias && os != NULL)
            {
                mp_real_set(os, x);
                mp_real_sin_cos_pi_bits(os, oc, os, prec);
            }
            else if (alias)
            {
                mp_real_set(oc, x);
                mp_real_sin_cos_pi_bits(os, oc, oc, prec);
            }
            else
                mp_real_sin_cos_pi_bits(os, oc, x, prec);

            if (os != NULL)
            {
                mp_real_get_arb(ya, os);
                trig_pi_test_check(ya, rs, xa, inexact, prec, "sin");
            }
            if (oc != NULL)
            {
                mp_real_get_arb(ya, oc);
                trig_pi_test_check(ya, rc, xa, inexact, prec, "cos");
            }
        }
        else
        {
            int ok;

            if (alias)
            {
                mp_real_set(y1, x);
                ok = mp_real_tan_pi_bits(y1, y1, prec);
            }
            else
                ok = mp_real_tan_pi_bits(y1, x, prec);

            if (ok)
            {
                mp_real_get_arb(ya, y1);
                trig_pi_test_check(ya, rt, xa, inexact, prec, "tan");
            }
            else
            {
                /* a pole in [m +- 3r] */
                arb_t w;
                arb_init(w);
                arb_set(w, xa);
                mag_mul_ui(arb_radref(w), arb_radref(w), 3);
                arb_cos_pi(t, w, wp);
                if (!arb_contains_zero(t) || !mp_real_is_zero(y1))
                {
                    flint_printf("FAIL: tan status 0, prec %wd\n", prec);
                    flint_printf("x = "); arb_printd(xa, 30); flint_printf("\n");
                    flint_abort();
                }
                arb_clear(w);
            }
        }

        mp_real_clear(x);
        mp_real_clear(y1);
        mp_real_clear(y2);
        arb_clear(xa);
        arb_clear(ya);
        arb_clear(rs);
        arb_clear(rc);
        arb_clear(rt);
        arb_clear(t);
        arf_clear(m);
        arf_clear(pt);
    }

    TEST_FUNCTION_END(state);
}
