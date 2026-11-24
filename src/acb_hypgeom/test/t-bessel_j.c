/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_bessel_j, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu0, nu1, nu2, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;

        acb_init(nu0);
        acb_init(nu1);
        acb_init(nu2);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_randtest_param(nu0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        acb_sub_ui(nu1, nu0, 1, prec0);
        acb_sub_ui(nu2, nu0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_j_asymp(w0, nu0, z, prec0);
                break;
            case 1:
                acb_hypgeom_bessel_j_0f1(w0, nu0, z, prec0);
                break;
            default:
                acb_hypgeom_bessel_j(w0, nu0, z, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_j_asymp(w1, nu0, z, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_j_0f1(w1, nu0, z, prec1);
                break;
            default:
                acb_hypgeom_bessel_j(w1, nu0, z, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("nu = "); acb_printd(nu0, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_j_asymp(w1, nu1, z, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_j_0f1(w1, nu1, z, prec1);
                break;
            default:
                acb_hypgeom_bessel_j(w1, nu1, z, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_j_asymp(w2, nu2, z, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_j_0f1(w2, nu2, z, prec2);
                break;
            default:
                acb_hypgeom_bessel_j(w2, nu2, z, prec2);
        }

        acb_mul(t, w1, nu1, prec0);
        acb_mul_2exp_si(t, t, 1);
        acb_submul(t, w2, z, prec0);
        acb_submul(t, w0, z, prec0);

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("nu = "); acb_printd(nu0, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_neg(t, nu0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_j_asymp(w2, t, z, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_j_0f1(w2, t, z, prec2);
                break;
            default:
                acb_hypgeom_bessel_j(w2, t, z, prec2);
        }

        acb_mul(w1, w1, w2, prec2);
        acb_neg(t, nu1);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_j_asymp(w2, t, z, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_j_0f1(w2, t, z, prec2);
                break;
            default:
                acb_hypgeom_bessel_j(w2, t, z, prec2);
        }

        acb_mul(w0, w0, w2, prec2);
        acb_add(w0, w0, w1, prec2);

        acb_sin_pi(t, nu0, prec2);
        acb_const_pi(u, prec2);
        acb_mul(u, u, z, prec2);
        acb_div(t, t, u, prec2);
        acb_mul_2exp_si(t, t, 1);

        if (!acb_overlaps(w0, t))
        {
            flint_printf("FAIL: wronskian\n\n");
            flint_printf("nu = "); acb_printd(nu0, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(nu0);
        acb_clear(nu1);
        acb_clear(nu2);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
        acb_clear(u);
    }

    /* Half-integer orders (the terminating asymptotic expansions where
       they are the cheapest method): J, I and K against the series at
       high precision, with nearly the full accuracy. */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, z, w, r;
        slong prec, rprec;
        int kind;

        acb_init(nu);
        acb_init(z);
        acb_init(w);
        acb_init(r);

        prec = 30 + n_randint(state, 400);
        rprec = 2 * prec + 1200;
        kind = n_randint(state, 3);

        acb_set_si(nu, 2 * ((slong) n_randint(state, 80) - 40) + 1);
        acb_mul_2exp_si(nu, nu, -1);

        /* |z| from 2^-10 to 2^8 */
        arb_set_si(acb_realref(z), 1 + n_randint(state, 1000));
        if (n_randint(state, 2))
            arb_set_si(acb_imagref(z), (slong) n_randint(state, 2000) - 1000);
        acb_mul_2exp_si(z, z, (slong) n_randint(state, 18) - 20);

        if (kind == 0)
        {
            acb_hypgeom_bessel_j(w, nu, z, prec);
            acb_hypgeom_bessel_j_0f1(r, nu, z, rprec);
        }
        else if (kind == 1)
        {
            acb_hypgeom_bessel_i(w, nu, z, prec);
            acb_hypgeom_bessel_i_0f1(r, nu, z, 0, rprec);
        }
        else
        {
            acb_hypgeom_bessel_k(w, nu, z, prec);
            acb_hypgeom_bessel_k_0f1(r, nu, z, 0, rprec);
        }

        if (!acb_overlaps(w, r) || acb_rel_accuracy_bits(r) < prec + 20
            || acb_rel_accuracy_bits(w) < prec - 16)
        {
            flint_printf("FAIL: half-integer order, kind %d\n\n", kind);
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(nu);
        acb_clear(z);
        acb_clear(w);
        acb_clear(r);
    }

    /* Inexact input: the result contains the values at points of the
       input balls, for J, Y, I and K with an exact or inexact order
       (also negative, where the bound for the variation with the
       order comes from the recurrences) and with an exact argument or
       a ball in the right half-plane; and close to a zero of J_0, J_1
       or Y_0 the radius is within some bits of that of the image of
       the ball. */
    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        static const char * zeros[7] = {
            "2.404825557695772768621631879", "5.520078110286310649596604113",
            "8.653727912911012216954198712", "3.831705970207512315614435886",
            "7.015586669815618753537049934", "0.893576966279167521584887102",
            "3.957678419314857868375677406" };
        static const int zero_kind[7] = { 0, 0, 0, 0, 0, 1, 1 };
        static const int zero_nu[7] = { 0, 0, 0, 1, 1, 0, 0 };
        acb_t nu, z, w, r, nup, zp;
        arf_t x;
        mag_t m;
        slong prec, rprec;
        int kind, k, zero = -1, mode;

        acb_init(nu);
        acb_init(z);
        acb_init(w);
        acb_init(r);
        acb_init(nup);
        acb_init(zp);
        arf_init(x);
        mag_init(m);

        prec = 30 + n_randint(state, 300);
        rprec = prec + 100;
        kind = n_randint(state, 4);
        mode = n_randint(state, 8);

        if (mode == 0)
        {
            zero = n_randint(state, 7);
            kind = zero_kind[zero];
            acb_set_si(nu, zero_nu[zero]);
            arb_set_str(acb_realref(z), zeros[zero], 90);
            acb_get_mid(z, z);
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -(slong) (10 + n_randint(state, prec)));
        }
        else if (mode <= 2)
        {
            /* Y and K at an integer order on a real ball (the
               derivative from the Wronskians) */
            kind = 1 + 2 * n_randint(state, 2);
            acb_set_si(nu, n_randint(state, 20));
            arb_set_si(acb_realref(z), 3 * 64 + n_randint(state, 40 * 64));
            acb_mul_2exp_si(z, z, -6);
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -(slong) (8 + n_randint(state, prec + 10)));
        }
        else if (mode == 3)
        {
            /* J and I with a negative order which is a ball (the
               bound from the recurrences) */
            kind = 2 * n_randint(state, 2);
            acb_set_si(nu, -1 - (slong) n_randint(state, 100));
            acb_mul_2exp_si(nu, nu, -3);
            mag_set_ui_2exp_si(arb_radref(acb_realref(nu)), 1, -(slong) (8 + n_randint(state, prec + 10)));
            arb_set_si(acb_realref(z), 1 + n_randint(state, 1024));
            acb_mul_2exp_si(z, z, -4);
        }
        else
        {
            /* the order: an integer or a multiple of 1/8, exact or a ball */
            acb_set_si(nu, (slong) n_randint(state, 300) - 80);
            if (n_randint(state, 2))
                acb_mul_2exp_si(nu, nu, -3);
            if (n_randint(state, 2))
                mag_set_ui_2exp_si(arb_radref(acb_realref(nu)), 1, -(slong) (8 + n_randint(state, prec + 10)));
            if (n_randint(state, 8) == 0)
                mag_set_ui_2exp_si(arb_radref(acb_imagref(nu)), 1, -(slong) (8 + n_randint(state, prec + 10)));

            /* the argument: 0 < Re(z) < 64, exact or a ball */
            arb_set_si(acb_realref(z), 1 + n_randint(state, 1024));
            if (n_randint(state, 3) == 0)
                arb_set_si(acb_imagref(z), (slong) n_randint(state, 1024) - 512);
            acb_mul_2exp_si(z, z, -4);
            if (n_randint(state, 3))
                mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -(slong) (8 + n_randint(state, prec + 10)));
            if (!arb_is_zero(acb_imagref(z)) && n_randint(state, 2))
                mag_set_ui_2exp_si(arb_radref(acb_imagref(z)), 1, -(slong) (8 + n_randint(state, prec + 10)));
        }

        switch (kind)
        {
            case 0: acb_hypgeom_bessel_j(w, nu, z, prec); break;
            case 1: acb_hypgeom_bessel_y(w, nu, z, prec); break;
            case 2: acb_hypgeom_bessel_i(w, nu, z, prec); break;
            default: acb_hypgeom_bessel_k(w, nu, z, prec);
        }

        for (k = 0; k < 3; k++)
        {
            /* a point of each ball: the midpoint plus the radius times
               -1, -1/2, 0, 1/2 or 1 */
            acb_get_mid(nup, nu);
            acb_get_mid(zp, z);

#define SHIFT(mid, rad) \
    do { \
        arf_set_mag(x, rad); \
        arf_mul_si(x, x, (slong) n_randint(state, 5) - 2, ARF_PREC_EXACT, ARF_RND_DOWN); \
        arf_mul_2exp_si(x, x, -1); \
        arf_add(mid, mid, x, ARF_PREC_EXACT, ARF_RND_DOWN); \
    } while (0)

            SHIFT(arb_midref(acb_realref(nup)), arb_radref(acb_realref(nu)));
            SHIFT(arb_midref(acb_imagref(nup)), arb_radref(acb_imagref(nu)));
            SHIFT(arb_midref(acb_realref(zp)), arb_radref(acb_realref(z)));
            SHIFT(arb_midref(acb_imagref(zp)), arb_radref(acb_imagref(z)));

            switch (kind)
            {
                case 0: acb_hypgeom_bessel_j(r, nup, zp, rprec); break;
                case 1: acb_hypgeom_bessel_y(r, nup, zp, rprec); break;
                case 2: acb_hypgeom_bessel_i(r, nup, zp, rprec); break;
                default: acb_hypgeom_bessel_k(r, nup, zp, rprec);
            }

            if (!acb_overlaps(w, r))
            {
                flint_printf("FAIL: point of the input balls, kind %d\n\n", kind);
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n\n");
                flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("nup = "); acb_printd(nup, 30); flint_printf("\n\n");
                flint_printf("zp = ");  acb_printd(zp, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (zero >= 0)
        {
            /* the derivative at the zero: -J_1, J_0 (as J_1 vanishes), -Y_1 */
            acb_get_mid(zp, z);
            if (zero_kind[zero] == 1)
            {
                acb_one(nup);
                acb_hypgeom_bessel_y(r, nup, zp, 64);
            }
            else
            {
                acb_set_si(nup, 1 - zero_nu[zero]);
                acb_hypgeom_bessel_j(r, nup, zp, 64);
            }

            /* rad(w) <= 2^7 (r |g'| + 2^-prec) */
            acb_get_mag(m, r);
            mag_mul(m, m, arb_radref(acb_realref(z)));
            mag_add_ui_2exp_si(m, m, 1, -prec);
            mag_mul_2exp_si(m, m, 7);

            if (mag_cmp(arb_radref(acb_realref(w)), m) > 0)
            {
                flint_printf("FAIL: radius close to a zero (%d)\n\n", zero);
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("bound = "); mag_printd(m, 10); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(nu);
        acb_clear(z);
        acb_clear(w);
        acb_clear(r);
        acb_clear(nup);
        acb_clear(zp);
        arf_clear(x);
        mag_clear(m);
    }

    /* A ball z which meets the cut (-inf, 0] from both sides: the
       values just above and just below the cut are in the result. */
    for (iter = 0; iter < 400 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, z, w, r, zp;
        slong prec;
        int kind, side;

        acb_init(nu);
        acb_init(z);
        acb_init(w);
        acb_init(r);
        acb_init(zp);

        prec = 2 + n_randint(state, 400);
        kind = n_randint(state, 4);

        /* nu = n/3 (an integer one time in three), or a ball */
        acb_set_si(nu, (slong) n_randint(state, 31) - 15);
        acb_div_ui(nu, nu, 3, prec);
        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(nu)), 1, -10 - (slong) n_randint(state, 40));

        /* z = -x + (+/- r) (1 + i), with 0 < x <= 30 and 2^-60 <= r <= 2^-4 */
        arb_set_ui(acb_realref(z), 1 + n_randint(state, 120));
        arb_mul_2exp_si(acb_realref(z), acb_realref(z), -2);
        arb_neg(acb_realref(z), acb_realref(z));
        arb_zero(acb_imagref(z));
        mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -4 - (slong) n_randint(state, 57));
        mag_set(arb_radref(acb_imagref(z)), arb_radref(acb_realref(z)));

        switch (kind)
        {
            case 0: acb_hypgeom_bessel_j(w, nu, z, prec); break;
            case 1: acb_hypgeom_bessel_y(w, nu, z, prec); break;
            case 2: acb_hypgeom_bessel_i(w, nu, z, prec); break;
            default: acb_hypgeom_bessel_k(w, nu, z, prec); break;
        }

        for (side = -1; side <= 1; side += 2)
        {
            /* the midpoint, moved by half the radius above or below */
            acb_get_mid(zp, z);
            arf_set_mag(arb_midref(acb_imagref(zp)), arb_radref(acb_imagref(z)));
            arf_mul_2exp_si(arb_midref(acb_imagref(zp)), arb_midref(acb_imagref(zp)), -1);
            if (side < 0)
                arf_neg(arb_midref(acb_imagref(zp)), arb_midref(acb_imagref(zp)));
            acb_get_mid(r, nu);

            switch (kind)
            {
                case 0: acb_hypgeom_bessel_j(r, r, zp, prec + 30); break;
                case 1: acb_hypgeom_bessel_y(r, r, zp, prec + 30); break;
                case 2: acb_hypgeom_bessel_i(r, r, zp, prec + 30); break;
                default: acb_hypgeom_bessel_k(r, r, zp, prec + 30); break;
            }

            if (!acb_overlaps(w, r))
            {
                flint_printf("FAIL: ball across the cut, kind %d\n\n", kind);
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n\n");
                flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("zp = ");  acb_printd(zp, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(nu);
        acb_clear(z);
        acb_clear(w);
        acb_clear(r);
        acb_clear(zp);
    }

    TEST_FUNCTION_END(state);
}
