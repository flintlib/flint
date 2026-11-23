/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "acb_hypgeom/sd.h"

/* The steepest-descent integration for large parameters, where the
   contours are planned and integrated in frames (sd.h):

   1. with frames forced for parameters of moderate size, against the
      integration in double precision coordinates;

   2. for parameters between 2^44 and 2^850, by identities between
      functions computed along different contours: the Wronskians of
      (I, K) and of the Hankel functions, and the recurrence of the
      incomplete gamma functions. */

static double
frnd(flint_rand_t state)
{
    return n_randint(state, 20001) / 10000.0 - 1.0;
}

static int
frames_integrate(acb_t res, const acb_t alpha, int sym, const acb_t nu, const acb_t start,
    sd_end a, sd_end b, slong prec, int frames)
{
    sd_problem P;
    sd_work W;
    sd_contour C;
    double tol = prec * 0.6931471805599453 + 12.0;
    int ok = 0;

    if (!_sd_problem_init_frames(&P, alpha, sym, nu, start, tol, tol, 16.0, frames))
        return 0;

    _sd_work_init(&W, &P);
    if (_sd_plan(&C, &W, &P, a, b, 0))
    {
        ok = _sd_integrate(res, &C, &P, alpha, sym, nu, start, prec);
        _sd_contour_clear(&C);
    }
    _sd_work_clear(&W);
    return ok;
}

/* a number of size about 2^e */
static void
rand_large(acb_t x, flint_rand_t state, slong e, int real)
{
    arb_set_d(acb_realref(x), (n_randint(state, 2) ? 1 : -1) * (0.5 + 0.5 * fabs(frnd(state))));
    arb_set_d(acb_imagref(x), real ? 0.0 : frnd(state));
    if (!real && n_randint(state, 4) == 0)
        arb_mul_2exp_si(acb_imagref(x), acb_imagref(x), -(slong) n_randint(state, 60));
    acb_mul_2exp_si(x, x, e);
}

TEST_FUNCTION_START(acb_hypgeom_sd_frames, state)
{
    slong iter, nboth = 0, nlarge = 0;

    /* frames against doubles */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        acb_t alpha, nu, start, r0, r1;
        sd_end a, b;
        slong prec = 30 + n_randint(state, 300);
        double s;
        int sym = n_randint(state, 2), ok0, ok1;

        acb_init(alpha); acb_init(nu); acb_init(start); acb_init(r0); acb_init(r1);

        s = exp(log(1e5) * n_randint(state, 1000) / 1000.0);
        acb_set_d_d(alpha, s * frnd(state), n_randint(state, 3) ? s * frnd(state) : 0.0);
        s = exp(log(1e5) * n_randint(state, 1000) / 1000.0);
        acb_set_d_d(nu, s * frnd(state), n_randint(state, 3) ? s * frnd(state) : 0.0);
        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(nu)), 1, -30 - (slong) n_randint(state, 100));
        if (sym && n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_imagref(alpha)), 1, -30 - (slong) n_randint(state, 100));

        if (sym)
        {
            a = sd_end_make(n_randint(state, 2) ? SD_PLUS : SD_MINUS, (int) n_randint(state, 3) - 1);
            do {
                b = sd_end_make(n_randint(state, 2) ? SD_PLUS : SD_MINUS, (int) n_randint(state, 3) - 1);
            } while (a.type == b.type && a.k == b.k);

            ok0 = frames_integrate(r0, alpha, 1, nu, NULL, a, b, prec, 0);
            ok1 = frames_integrate(r1, alpha, 1, nu, NULL, a, b, prec, 1);
        }
        else
        {
            acb_set_si(alpha, -1);
            acb_set_d_d(start, log(s) * (n_randint(state, 1000) / 500.0), 2.5 * frnd(state));
            if (n_randint(state, 3) == 0)
                mag_set_ui_2exp_si(arb_radref(acb_realref(start)), 1, -30 - (slong) n_randint(state, 100));
            a = sd_end_make(SD_START, 0);
            b = sd_end_make(n_randint(state, 2) ? SD_PLUS : SD_LINEAR, 0);

            ok0 = frames_integrate(r0, alpha, 0, nu, start, a, b, prec, 0);
            ok1 = frames_integrate(r1, alpha, 0, nu, start, a, b, prec, 1);
        }

        nboth += ok0 && ok1;

        /* the same contours are found, and the results agree; with
           frames the accuracy is as good, or the precision is reached */
        if ((ok0 && !ok1) || (ok0 && ok1 && (!acb_overlaps(r0, r1)
            || (acb_rel_accuracy_bits(r1) < acb_rel_accuracy_bits(r0) - 10
                && acb_rel_accuracy_bits(r1) < prec - 4))))
        {
            flint_printf("FAIL: frames against doubles\n\n");
            flint_printf("sym = %d, prec = %wd, ok = %d, %d, ends (%d, %d), (%d, %d)\n",
                sym, prec, ok0, ok1, a.type, a.k, b.type, b.k);
            flint_printf("alpha = "); acb_printd(alpha, 30); flint_printf("\n");
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
            flint_printf("start = "); acb_printd(start, 30); flint_printf("\n");
            flint_printf("r0 = "); acb_printd(r0, 30); flint_printf("\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(alpha); acb_clear(nu); acb_clear(start); acb_clear(r0); acb_clear(r1);
    }

    /* large parameters */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        acb_t nu, nu1, z, a0, a1, b0, b1, t, u;
        slong prec = 40 + n_randint(state, 200), e = 44 + n_randint(state, 800), acc, wp;
        int real = n_randint(state, 2), type = n_randint(state, 4), which = n_randint(state, 3);
        int exact = n_randint(state, 4) != 0, upper = n_randint(state, 2), ok = 1;

        acb_init(nu); acb_init(nu1); acb_init(z); acb_init(a0); acb_init(a1);
        acb_init(b0); acb_init(b1); acb_init(t); acb_init(u);

        /* nu of the size of z, close to the turning point nu = z (or
           nu = i z), smaller, or larger */
        rand_large(z, state, e, real);
        if (type == 0)
            rand_large(nu, state, e + (slong) n_randint(state, 5) - 2, real);
        else if (type == 1)
        {
            acb_set(nu, z);
            if (!real && n_randint(state, 2))
                acb_mul_onei(nu, nu);
            arb_set_d(acb_realref(t), 1.0 + 1000.0 * frnd(state) * ldexp(1.0, -(int) n_randint(state, 50)));
            arb_zero(acb_imagref(t));
            acb_mul(nu, nu, t, 2000);
        }
        else if (type == 2)
            rand_large(nu, state, e - (slong) n_randint(state, 60), real);
        else
            rand_large(nu, state, FLINT_MIN(e + (slong) n_randint(state, 30), 850), real);

        if (real && arb_is_negative(acb_realref(z)))
            acb_neg(z, z);

        acb_add_ui(nu1, nu, 1, 4000);

        if (!exact)
        {
            acb_get_mag(arb_radref(acb_realref(t)), z);
            mag_mul_2exp_si(arb_radref(acb_realref(z)), arb_radref(acb_realref(t)), -prec / 2 - 20);
        }

        if (which == 0)
        {
            /* z (I_nu K_{nu+1} + I_{nu+1} K_nu) = 1 */
            acb_hypgeom_bessel_i_integration(a0, nu, z, prec);
            acb_hypgeom_bessel_i_integration(a1, nu1, z, prec);
            acb_hypgeom_bessel_k_integration(b0, nu, z, prec);
            acb_hypgeom_bessel_k_integration(b1, nu1, z, prec);
            acb_mul(t, a0, b1, prec);
            acb_addmul(t, a1, b0, prec);
            acb_mul(t, t, z, prec);
            acb_sub_ui(t, t, 1, prec);
        }
        else if (which == 1)
        {
            /* H1_{nu+1} H2_nu - H1_nu H2_{nu+1} = -4 i / (pi z) */
            acb_hypgeom_bessel_h1_integration(a0, nu, z, prec);
            acb_hypgeom_bessel_h1_integration(a1, nu1, z, prec);
            acb_hypgeom_bessel_h2_integration(b0, nu, z, prec);
            acb_hypgeom_bessel_h2_integration(b1, nu1, z, prec);
            acb_mul(t, a1, b0, prec);
            acb_submul(t, a0, b1, prec);
            acb_mul(t, t, z, prec);
            acb_const_pi(u, prec);
            acb_mul(t, t, u, prec);
            acb_mul_2exp_si(t, t, -2);
            acb_div_onei(t, t);
            acb_add_ui(t, t, 1, prec);
        }
        else
        {
            /* Gamma(s+1, z) - s Gamma(s, z) = z^s e^-z, and the lower
               function with the opposite sign (which does not exist
               for the negative integers s among the large doubles),
               for z below 2^850 */
            if (e > 840)
                acb_mul_2exp_si(z, z, -30);
            ok = _acb_hypgeom_gamma_incomplete_integration(a0, nu, z, upper, prec)
                && _acb_hypgeom_gamma_incomplete_integration(a1, nu1, z, upper, prec);

            if (ok)
            {
                wp = prec + 2 * e + 100;
                acb_log(u, z, wp);
                acb_mul(u, u, nu, wp);
                acb_sub(u, u, z, wp);
                acb_neg(u, u);
                acb_exp(u, u, wp);
                acb_mul(t, nu, a0, wp);
                acb_sub(t, a1, t, wp);
                acb_mul(t, t, u, wp);
                if (upper)
                    acb_sub_ui(t, t, 1, wp);
                else
                    acb_add_ui(t, t, 1, wp);
                acb_set(b0, a0);
                acb_set(b1, a1);
            }
        }

        /* the integration applies (for a ball z, the functions vary
           too much to be enclosed) */
        if (!ok && exact && !(!upper && acb_is_int(nu) && arb_is_negative(acb_realref(nu))))
        {
            flint_printf("FAIL: large parameters (not integrated)\n\n");
            flint_printf("which = %d, upper = %d, prec = %wd\n", which, upper, prec);
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_abort();
        }

        if (ok)
        {
            nlarge++;

            acc = FLINT_MIN(acb_rel_accuracy_bits(a0), acb_rel_accuracy_bits(a1));
            acc = FLINT_MIN(acc, acb_rel_accuracy_bits(b0));
            acc = FLINT_MIN(acc, acb_rel_accuracy_bits(b1));

            if (!acb_contains_zero(t) || (exact && acc < prec - 16))
            {
                flint_printf("FAIL: large parameters\n\n");
                flint_printf("which = %d, upper = %d, prec = %wd, acc = %wd\n", which, upper, prec, acc);
                flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n");
                flint_printf("a1 = "); acb_printd(a1, 30); flint_printf("\n");
                flint_printf("b0 = "); acb_printd(b0, 30); flint_printf("\n");
                flint_printf("b1 = "); acb_printd(b1, 30); flint_printf("\n");
                flint_printf("t = "); acb_printd(t, 30); flint_printf("\n");
                flint_abort();
            }
        }

        acb_clear(nu); acb_clear(nu1); acb_clear(z); acb_clear(a0); acb_clear(a1);
        acb_clear(b0); acb_clear(b1); acb_clear(t); acb_clear(u);
    }

    /* (most of the cases must have been integrated) */
    if (nboth < 100 * flint_test_multiplier() || nlarge < 20 * flint_test_multiplier())
    {
        flint_printf("FAIL: %wd and %wd cases integrated\n", nboth, nlarge);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
