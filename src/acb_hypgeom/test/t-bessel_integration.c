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
#include "double_extras.h"
#include "acb_hypgeom/sd.h"

#define KIND_J ACB_HYPGEOM_KIND_BESSEL_J
#define KIND_Y ACB_HYPGEOM_KIND_BESSEL_Y
#define KIND_JY ACB_HYPGEOM_KIND_BESSEL_JY
#define KIND_H1 ACB_HYPGEOM_KIND_HANKEL_1
#define KIND_H2 ACB_HYPGEOM_KIND_HANKEL_2
#define KIND_I ACB_HYPGEOM_KIND_BESSEL_I
#define KIND_K ACB_HYPGEOM_KIND_BESSEL_K

static const int kinds[7] = { KIND_J, KIND_Y, KIND_JY, KIND_H1, KIND_H2, KIND_I, KIND_K };

/* a number of magnitude between 2^-6 and 2^maxexp: real, on a ray close
   to an axis, or complex; sometimes an integer or a half-integer */
static void
rand_exact(acb_t x, flint_rand_t state, int maxexp)
{
    double a, b = 0.0;

    a = ldexp(1.0 + n_randint(state, 1000) / 1000.0, (int) n_randint(state, maxexp + 6) - 6);
    if (n_randint(state, 2))
        a = -a;
    if (n_randint(state, 8) == 0)
        a = floor(a) + n_randint(state, 2) * 0.5;

    switch (n_randint(state, 6))
    {
        case 0:
            b = ldexp(1.0 + n_randint(state, 1000) / 1000.0, (int) n_randint(state, maxexp + 6) - 6);
            if (n_randint(state, 2))
                b = -b;
            break;
        case 1:
            b = a * ldexp(1.0, -(int) n_randint(state, 40));
            break;
        case 2:
            b = a;
            a = b * ldexp(1.0, -(int) n_randint(state, 40));
            break;
        default:
            break;
    }

    acb_set_d_d(x, a, b);
}

/* widen x to a ball (sometimes crossing the real axis), and set p to a
   random point of the ball */
static void
rand_ball(acb_t p, acb_t x, flint_rand_t state)
{
    arb_ptr c, q;
    arf_t t;
    mag_t r;
    int i;

    arf_init(t);
    mag_init(r);

    for (i = 0; i < 2; i++)
    {
        c = i ? acb_imagref(x) : acb_realref(x);
        q = i ? acb_imagref(p) : acb_realref(p);

        arb_set(q, c);

        if (n_randint(state, 2))
        {
            acb_get_mag(r, x);
            mag_mul_2exp_si(r, r, -2 - (slong) n_randint(state, 40));
            mag_set(arb_radref(c), r);
            arf_set_mag(t, r);
            arf_mul_si(t, t, (slong) n_randint(state, 201) - 100, 64, ARF_RND_DOWN);
            arf_div_ui(t, t, 101, 64, ARF_RND_DOWN);
            arb_add_arf(q, q, t, ARF_PREC_EXACT);
            mag_zero(arb_radref(q));
        }
    }

    arf_clear(t);
    mag_clear(r);
}

/* the reference values, by the series and asymptotic expansions */
static void
reference(acb_t res, acb_t res2, int kind, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    switch (kind)
    {
        case KIND_J:
            acb_hypgeom_bessel_j_nointegration(res, nu, z, prec);
            break;
        case KIND_Y:
            acb_hypgeom_bessel_y_nointegration(res, nu, z, prec);
            break;
        case KIND_I:
            _acb_hypgeom_bessel_i_nointegration(res, nu, z, scaled, prec);
            break;
        case KIND_K:
            _acb_hypgeom_bessel_k_nointegration(res, nu, z, scaled, prec);
            break;
        default:
            acb_hypgeom_bessel_jy_nointegration(res, res2, nu, z, prec);
            if (kind != KIND_JY)
            {
                acb_mul_onei(res2, res2);
                if (kind == KIND_H1)
                    acb_add(res, res, res2, prec);
                else
                    acb_sub(res, res, res2, prec);
            }
    }
}

static int
usable(const acb_t x)
{
    return acb_is_finite(x) && acb_rel_accuracy_bits(x) > 10;
}

static slong
reference_prec(const acb_t nu, const acb_t z, slong prec)
{
    mag_t t;
    slong p;

    mag_init(t);
    acb_get_mag(t, nu);
    p = 2 * prec + 100 + (slong) (4 * mag_get_d(t));
    acb_get_mag(t, z);
    p += (slong) (4 * mag_get_d(t));
    mag_clear(t);
    return p;
}

TEST_FUNCTION_START(acb_hypgeom_bessel_integration, state)
{
    slong iter;

    /* Enclosures: compare with the series and asymptotic expansions at
       a point of the input balls.  Mostly at the precisions where the
       quadrature uses doubles; sometimes higher. */
    for (iter = 0; iter < 400 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, z, nu0, z0, r, r2, s, s2;
        slong prec, rprec;
        int kind, scaled, ok;

        acb_init(nu); acb_init(z); acb_init(nu0); acb_init(z0);
        acb_init(r); acb_init(r2); acb_init(s); acb_init(s2);

        kind = kinds[n_randint(state, 7)];
        scaled = (kind == KIND_I || kind == KIND_K) && n_randint(state, 2);
        prec = 10 + n_randint(state, (n_randint(state, 10) == 0) ? 500 : 190);

        rand_exact(nu, state, 5);
        rand_exact(z, state, 5);
        acb_set(nu0, nu);
        acb_set(z0, z);

        if (iter % 2)
        {
            rand_ball(nu0, nu, state);
            rand_ball(z0, z, state);
        }

        rprec = reference_prec(nu, z, prec);
        reference(s, s2, kind, nu0, z0, scaled, rprec);

        /* also with aliasing */
        if (n_randint(state, 4) == 0)
        {
            acb_set(r, nu);
            ok = _acb_hypgeom_bessel_integration(r, r2, kind, r, z, scaled, prec);
        }
        else
            ok = _acb_hypgeom_bessel_integration(r, r2, kind, nu, z, scaled, prec);

        if (ok && ((usable(s) && !acb_overlaps(r, s))
            || (kind == KIND_JY && usable(s2) && !acb_overlaps(r2, s2))))
        {
            flint_printf("FAIL: enclosure\n\n");
            flint_printf("kind = %d, scaled = %d, prec = %wd\n", kind, scaled, prec);
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("nu0 = "); acb_printd(nu0, 30); flint_printf("\n");
            flint_printf("z0 = "); acb_printd(z0, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n");
            flint_printf("s2 = "); acb_printd(s2, 30); flint_printf("\n");
            flint_abort();
        }

        /* exact, real and positive: the contours must be found, and
           there is no cancellation */
        if (acb_is_exact(nu) && acb_is_exact(z) && acb_is_real(nu) && acb_is_real(z)
            && arb_is_positive(acb_realref(nu)) && arb_is_positive(acb_realref(z))
            && usable(s) && (!ok || acb_rel_accuracy_bits(r) < prec - 12))
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("kind = %d, scaled = %d, prec = %wd, ok = %d\n", kind, scaled, prec, ok);
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(nu); acb_clear(z); acb_clear(nu0); acb_clear(z0);
        acb_clear(r); acb_clear(r2); acb_clear(s); acb_clear(s2);
    }

    /* Large parameters, where the series are of no use: the Wronskians
       J_{nu+1} Y_nu - J_nu Y_{nu+1} = 2 / (pi z),
       I_nu K_{nu+1} + I_{nu+1} K_nu = 1 / z, and H1,2 = J +- i Y */
    for (iter = 0; iter < 40 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, nu1, z, a, b, c, d, t, u;
        slong prec;
        int ok, bessel_ik;

        acb_init(nu); acb_init(nu1); acb_init(z);
        acb_init(a); acb_init(b); acb_init(c); acb_init(d);
        acb_init(t); acb_init(u);

        prec = 20 + n_randint(state, 150);
        rand_exact(nu, state, 2 + n_randint(state, 30));
        rand_exact(z, state, 2 + n_randint(state, 30));
        acb_add_ui(nu1, nu, 1, ARF_PREC_EXACT);
        bessel_ik = n_randint(state, 2);

        if (bessel_ik)
        {
            ok = _acb_hypgeom_bessel_integration(a, NULL, KIND_I, nu, z, 0, prec)
                && _acb_hypgeom_bessel_integration(b, NULL, KIND_K, nu1, z, 0, prec)
                && _acb_hypgeom_bessel_integration(c, NULL, KIND_I, nu1, z, 0, prec)
                && _acb_hypgeom_bessel_integration(d, NULL, KIND_K, nu, z, 0, prec);
            acb_mul(t, a, b, prec);
            acb_addmul(t, c, d, prec);
            acb_inv(u, z, prec);
        }
        else
        {
            /* (J, Y) = (a, b) and (c, d) */
            ok = _acb_hypgeom_bessel_integration(a, b, KIND_JY, nu, z, 0, prec)
                && _acb_hypgeom_bessel_integration(c, d, KIND_JY, nu1, z, 0, prec);
            acb_mul(t, c, b, prec);
            acb_submul(t, a, d, prec);
            acb_const_pi(u, prec);
            acb_mul(u, u, z, prec);
            acb_mul_2exp_si(u, u, -1);
            acb_inv(u, u, prec);
        }

        if (ok && !acb_overlaps(t, u))
        {
            flint_printf("FAIL: Wronskian\n\n");
            flint_printf("ik = %d, prec = %wd\n", bessel_ik, prec);
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n");
            flint_printf("u = "); acb_printd(u, 30); flint_printf("\n");
            flint_abort();
        }

        if (ok && !bessel_ik)
        {
            int kind = n_randint(state, 2) ? KIND_H1 : KIND_H2;

            acb_mul_onei(t, b);
            if (kind == KIND_H1)
                acb_add(t, a, t, prec);
            else
                acb_sub(t, a, t, prec);

            if (_acb_hypgeom_bessel_integration(u, NULL, kind, nu, z, 0, prec)
                && !acb_overlaps(t, u))
            {
                flint_printf("FAIL: Hankel\n\n");
                flint_printf("kind = %d, prec = %wd\n", kind, prec);
                flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                flint_printf("t = "); acb_printd(t, 30); flint_printf("\n");
                flint_printf("u = "); acb_printd(u, 30); flint_printf("\n");
                flint_abort();
            }
        }

        acb_clear(nu); acb_clear(nu1); acb_clear(z);
        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(d);
        acb_clear(t); acb_clear(u);
    }

    /* The bounds: |F| <= M and |F(nu1, z1) - F(nu0, z0)| <= Dz |z1 - z0|
       + Dnu |nu1 - nu0| at two points of the input balls */
    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, z, nu0, z0, nu1, z1, f0, f1, t;
        mag_t M, Dz, Dnu, lhs, rhs, m;
        slong prec = 100;
        int kind;

        acb_init(nu); acb_init(z); acb_init(nu0); acb_init(z0);
        acb_init(nu1); acb_init(z1); acb_init(f0); acb_init(f1); acb_init(t);
        mag_init(M); mag_init(Dz); mag_init(Dnu); mag_init(lhs); mag_init(rhs); mag_init(m);

        kind = kinds[n_randint(state, 7)];
        if (kind == KIND_JY || kind == KIND_H1 || kind == KIND_H2)
            kind = KIND_J;

        rand_exact(nu, state, 5);
        rand_exact(z, state, 5);
        acb_set(nu0, nu);
        acb_set(z0, z);
        rand_ball(nu1, nu, state);
        rand_ball(z1, z, state);

        if (_acb_hypgeom_bessel_integration_bounds(M, Dz, Dnu, kind, nu, z))
        {
            reference(f0, t, kind, nu0, z0, 0, reference_prec(nu, z, prec));
            reference(f1, t, kind, nu1, z1, 0, reference_prec(nu, z, prec));

            if (usable(f0) && usable(f1))
            {
                acb_get_mag_lower(lhs, f0);

                if (mag_cmp(lhs, M) > 0)
                {
                    flint_printf("FAIL: magnitude bound\n\n");
                    flint_printf("kind = %d\n", kind);
                    flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                    flint_printf("f0 = "); acb_printd(f0, 30); flint_printf("\n");
                    flint_printf("M = "); mag_printd(M, 10); flint_printf("\n");
                    flint_abort();
                }

                acb_sub(t, f1, f0, prec);
                acb_get_mag_lower(lhs, t);
                acb_sub(t, z1, z0, prec);
                acb_get_mag(m, t);
                mag_mul(rhs, Dz, m);
                acb_sub(t, nu1, nu0, prec);
                acb_get_mag(m, t);
                mag_addmul(rhs, Dnu, m);

                if (mag_cmp(lhs, rhs) > 0)
                {
                    flint_printf("FAIL: derivative bounds\n\n");
                    flint_printf("kind = %d\n", kind);
                    flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                    flint_printf("nu1 = "); acb_printd(nu1, 30); flint_printf("\n");
                    flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
                    flint_printf("Dz = "); mag_printd(Dz, 10); flint_printf("\n");
                    flint_printf("Dnu = "); mag_printd(Dnu, 10); flint_printf("\n");
                    flint_printf("lhs = "); mag_printd(lhs, 10); flint_printf("\n");
                    flint_printf("rhs = "); mag_printd(rhs, 10); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        acb_clear(nu); acb_clear(z); acb_clear(nu0); acb_clear(z0);
        acb_clear(nu1); acb_clear(z1); acb_clear(f0); acb_clear(f1); acb_clear(t);
        mag_clear(M); mag_clear(Dz); mag_clear(Dnu); mag_clear(lhs); mag_clear(rhs); mag_clear(m);
    }

    /* The automatic functions with inexact input (evaluation at the
       midpoint with derivative bounds, among other routes): compare
       with the value at a point of the input balls; often with only z
       inexact, and with balls crossing the branch cut */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu, z, nu0, z0, r, f, f2;
        slong prec;
        int kind, scaled;

        acb_init(nu); acb_init(z); acb_init(nu0); acb_init(z0);
        acb_init(r); acb_init(f); acb_init(f2);

        kind = kinds[n_randint(state, 2) ? n_randint(state, 2) : 5 + n_randint(state, 2)];
        scaled = (kind == KIND_I || kind == KIND_K) && n_randint(state, 2);
        prec = 10 + n_randint(state, 190);
        rand_exact(nu, state, 5);
        rand_exact(z, state, 5);

        if (n_randint(state, 3) == 0)
            acb_set_si(nu, (slong) n_randint(state, 20) - 5);

        if (n_randint(state, 4) == 0)
        {
            arb_zero(acb_imagref(z));
            if (n_randint(state, 2))
                arb_neg(acb_realref(z), acb_realref(z));
        }

        acb_set(nu0, nu);
        acb_set(z0, z);
        if (n_randint(state, 3) == 0)
            rand_ball(nu0, nu, state);
        rand_ball(z0, z, state);

        reference(f, f2, kind, nu0, z0, scaled, reference_prec(nu, z, prec));

        switch (kind)
        {
            case KIND_J: acb_hypgeom_bessel_j(r, nu, z, prec); break;
            case KIND_Y: acb_hypgeom_bessel_y(r, nu, z, prec); break;
            case KIND_I:
                if (scaled)
                    acb_hypgeom_bessel_i_scaled(r, nu, z, prec);
                else
                    acb_hypgeom_bessel_i(r, nu, z, prec);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(r, nu, z, prec);
                else
                    acb_hypgeom_bessel_k(r, nu, z, prec);
        }

        if (usable(f) && !acb_overlaps(r, f))
        {
            flint_printf("FAIL: enclosure (automatic)\n\n");
            flint_printf("kind = %d, scaled = %d, prec = %wd\n", kind, scaled, prec);
            flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("nu0 = "); acb_printd(nu0, 30); flint_printf("\n");
            flint_printf("z0 = "); acb_printd(z0, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
            flint_printf("f = "); acb_printd(f, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(nu); acb_clear(z); acb_clear(nu0); acb_clear(z0);
        acb_clear(r); acb_clear(f); acb_clear(f2);
    }

    /* The saddle-point expansions against the quadrature on the same
       contours, for the symmetric kernel with large parameters. */
    {
        slong nexp = 0, niter = 300 * 0.1 * flint_test_multiplier();

        for (iter = 0; iter < niter; iter++)
        {
            acb_t alpha, nu, r1, r2;
            sd_problem P;
            sd_work W;
            sd_contour C;
            sd_end a, b;
            slong prec;
            double tol, m, th;
            int ok1, ok2;

            acb_init(alpha); acb_init(nu); acb_init(r1); acb_init(r2);

            prec = 20 + n_randint(state, n_randint(state, 20) ? 150 : 1000);
            tol = prec * SD_LN2 + 12.0;

            m = ldexp(1.0 + d_randtest(state), 2 + (int) n_randint(state, 34));
            th = n_randint(state, 3) ? 6.3 * d_randtest(state) : 1.5707963267948966 * n_randint(state, 4);
            acb_set_d_d(alpha, m * cos(th), n_randint(state, 4) ? m * sin(th) : 0.0);
            m = ldexp(1.0 + d_randtest(state), 2 + (int) n_randint(state, 34));
            th = n_randint(state, 3) ? 6.3 * d_randtest(state) : 1.5707963267948966 * n_randint(state, 4);
            acb_set_d_d(nu, m * cos(th), n_randint(state, 4) ? m * sin(th) : 0.0);

            /* close to a double saddle (nu = +- 2 alpha) */
            if (n_randint(state, 4) == 0)
            {
                acb_set_d_d(r1, ldexp(d_randtest(state), -(int) n_randint(state, 20)),
                    n_randint(state, 2) ? ldexp(d_randtest(state), -(int) n_randint(state, 20)) : 0.0);
                acb_add_ui(r1, r1, 1, 53);
                acb_get_mid(r1, r1);
                acb_mul(nu, alpha, r1, 53);
                acb_get_mid(nu, nu);
                acb_mul_2exp_si(nu, nu, 1);
                if (n_randint(state, 2))
                    acb_neg(nu, nu);
            }

            a = sd_end_make(n_randint(state, 2) ? SD_PLUS : SD_MINUS, (int) n_randint(state, 3) - 1);
            do {
                b = sd_end_make(n_randint(state, 2) ? SD_PLUS : SD_MINUS, (int) n_randint(state, 3) - 1);
            } while (a.type == b.type && a.k == b.k);

            if (_sd_problem_init(&P, alpha, 1, nu, NULL, tol, tol, 16.0))
            {
                _sd_work_init(&W, &P);

                if (_sd_plan(&C, &W, &P, a, b, 1))
                {
                    ok1 = _sd_expand(r1, &C, &P, alpha, 1, nu, prec, D_INF);
                    ok2 = ok1 && _sd_integrate(r2, &C, &P, alpha, 1, nu, NULL, prec);
                    nexp += ok1;

                    if (ok1 && ok2 && !acb_overlaps(r1, r2))
                    {
                        flint_printf("FAIL: expansion\n\n");
                        flint_printf("prec = %wd, ends = (%d, %d), (%d, %d)\n", prec, a.type, a.k, b.type, b.k);
                        flint_printf("alpha = "); acb_printd(alpha, 30); flint_printf("\n");
                        flint_printf("nu = "); acb_printd(nu, 30); flint_printf("\n");
                        flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n");
                        flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n");
                        flint_abort();
                    }

                    _sd_contour_clear(&C);
                }

                _sd_work_clear(&W);
            }

            acb_clear(alpha); acb_clear(nu); acb_clear(r1); acb_clear(r2);
        }

        if (niter >= 300 && nexp < niter / 20)
        {
            flint_printf("FAIL: the expansion was used %wd times of %wd\n", nexp, niter);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
