/*
    Copyright (C) 2014 Fredrik Johansson

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

TEST_FUNCTION_START(acb_hypgeom_bessel_k, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu0, nu1, nu2, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;
        int scaled;

        acb_init(nu0);
        acb_init(nu1);
        acb_init(nu2);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 700);
        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        scaled = n_randint(state, 2);
        acb_randtest_param(nu0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(z));

        if (n_randint(state, 2) && prec0 <= 200 && prec1 <= 200 && prec2 <= 200)
        {
            arb_zero(acb_imagref(nu0));
            arb_zero(acb_imagref(z));
            arb_abs(acb_realref(nu0), acb_realref(nu0));
        }

        acb_sub_ui(nu1, nu0, 1, prec0);
        acb_sub_ui(nu2, nu0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w0, nu0, z, scaled, prec0);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w0, nu0, z, scaled, prec0);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(w0, nu0, z, prec0);
                else
                    acb_hypgeom_bessel_k(w0, nu0, z, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w1, nu0, z, scaled, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w1, nu0, z, scaled, prec1);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(w1, nu0, z, prec1);
                else
                    acb_hypgeom_bessel_k(w1, nu0, z, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("scaled = %d\n\n", scaled);
            flint_printf("nu = "); acb_printd(nu0, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w1, nu1, z, scaled, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w1, nu1, z, scaled, prec1);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(w1, nu1, z, prec1);
                else
                    acb_hypgeom_bessel_k(w1, nu1, z, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w2, nu2, z, scaled, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w2, nu2, z, scaled, prec2);
                break;
            default:
                if (scaled)
                {
                    acb_hypgeom_bessel_k(w2, nu2, z, prec2);
                    acb_exp(t, z, prec2);
                    acb_mul(w2, w2, t, prec2);
                }
                else
                {
                    acb_hypgeom_bessel_k_scaled(w2, nu2, z, prec2);
                    acb_neg(t, z);
                    acb_exp(t, t, prec2);
                    acb_mul(w2, w2, t, prec2);
                }
        }

        acb_mul(t, w1, nu1, prec0);
        acb_mul_2exp_si(t, t, 1);
        acb_addmul(t, w2, z, prec0);
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

    /* The sizes of I and K with complex orders, by which the
       cancellation of their series is estimated: |I_nu(z)| ~ e^eta and
       |K_nu(z)| ~ e^(-eta) up to algebraic factors, for Re(nu) >= 0,
       Re(z) > 0 and |arg(z) - arg(nu)| <= pi / 2.  (K is smaller close
       to its zeros, which lie on the real axis for imaginary orders.)
       Fixed cases: K_{i tau}(tau) is of the size e^(-pi tau / 2), where
       its series loses 0.36 tau bits, and is computed by the series. */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, z, w;
        arb_t m;
        double num, zn, th, ph, nr, ni, zr, zi, eta, lk, li, tol;

        acb_init(nu);
        acb_init(z);
        acb_init(w);
        arb_init(m);

        num = exp(log(5.0) + log(400.0) * n_randint(state, 1001) / 1000.0);
        zn = num * exp(log(0.02) + log(2500.0) * n_randint(state, 1001) / 1000.0);
        th = 1.5707963267948966 * (n_randint(state, 2001) / 1000.0 - 1.0);
        ph = 1.55 * (n_randint(state, 2001) / 1000.0 - 1.0);
        if (n_randint(state, 3) == 0)
            ph = 0.0;
        if (n_randint(state, 5) == 0)
            th = (n_randint(state, 2) ? 1 : -1) * 1.5707963267948966;
        if (n_randint(state, 6) == 0)
            zn = num;

        nr = floor(num * cos(th) * 16.0) / 16.0;
        nr = FLINT_MAX(nr, 0.0);
        ni = floor(num * sin(th) * 16.0) / 16.0;
        zr = floor(zn * cos(ph) * 16.0 + 1.0) / 16.0;
        zi = floor(zn * sin(ph) * 16.0) / 16.0;

        if (nr * zr + ni * zi >= 0.0)
        {
            acb_set_d_d(nu, nr, ni);
            acb_set_d_d(z, zr, zi);
            eta = _acb_hypgeom_bessel_debye_eta_re(nr, ni, zr, zi);
            tol = 5.0 + 0.25 * log(num * num + zn * zn);

            acb_hypgeom_bessel_k(w, nu, z, 128);
            acb_abs(m, w, 128);
            arb_log(m, m, 128);
            lk = arf_get_d(arb_midref(m), ARF_RND_NEAR);

            acb_hypgeom_bessel_i(w, nu, z, 128);
            acb_abs(m, w, 128);
            arb_log(m, m, 128);
            li = arf_get_d(arb_midref(m), ARF_RND_NEAR);

            if (!(fabs(li - eta) <= tol) || !(lk <= -eta + tol)
                || (!(nr == 0.0 && zi == 0.0) && !(lk >= -eta - tol)))
            {
                flint_printf("FAIL: sizes of I and K with a complex order\n\n");
                flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n\n");
                flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("log|I| = %g, log|K| = %g, eta = %g, tol = %g\n\n", li, lk, eta, tol);
                flint_abort();
            }
        }

        acb_clear(nu);
        acb_clear(z);
        acb_clear(w);
        arb_clear(m);
    }

    {
        acb_t nu, z, w, w2;
        slong i;

        acb_init(nu);
        acb_init(z);
        acb_init(w);
        acb_init(w2);

        for (i = 0; i < 3; i++)
        {
            slong prec = (i == 0) ? 333 : 2000;
            double tau = (i == 2) ? 3000.0 : 1000.0;

            acb_set_d_d(nu, 0.0, tau);
            acb_set_d(z, (i == 1) ? 500.0 : tau);

            acb_hypgeom_bessel_k(w, nu, z, prec);
            acb_hypgeom_bessel_k_0f1(w2, nu, z, 0, prec + (slong) (0.4 * tau) + 60);

            if (acb_rel_accuracy_bits(w) < prec - 10 || acb_rel_accuracy_bits(w2) < prec
                || !acb_overlaps(w, w2))
            {
                flint_printf("FAIL: K with an imaginary order, case %wd\n\n", i);
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(nu);
        acb_clear(z);
        acb_clear(w);
        acb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}
