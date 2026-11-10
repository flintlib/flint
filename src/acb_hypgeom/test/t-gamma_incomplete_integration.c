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

/* a number of magnitude between 2^-6 and 2^maxexp: real, on a ray close
   to an axis, or complex; sometimes an integer or a half-integer */
static void
gamma_rand_exact(acb_t x, flint_rand_t state, int maxexp)
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
gamma_rand_ball(acb_t p, acb_t x, flint_rand_t state)
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

static void
gamma_reference(acb_t res, const acb_t s, const acb_t z, int upper, slong prec)
{
    if (upper)
        acb_hypgeom_gamma_upper_nointegration(res, s, z, 0, prec);
    else
        acb_hypgeom_gamma_lower_nointegration(res, s, z, 0, prec);
}

static int
gamma_usable(const acb_t x)
{
    return acb_is_finite(x) && acb_rel_accuracy_bits(x) > 10;
}

static slong
gamma_reference_prec(const acb_t s, const acb_t z, slong prec)
{
    mag_t t;
    slong p;

    mag_init(t);
    acb_get_mag(t, s);
    p = 2 * prec + 100 + (slong) (4 * mag_get_d(t));
    acb_get_mag(t, z);
    p += (slong) (4 * mag_get_d(t));
    mag_clear(t);
    return p;
}

TEST_FUNCTION_START(acb_hypgeom_gamma_incomplete_integration, state)
{
    slong iter;

    /* Enclosures: compare with the series and asymptotic expansions at
       a point of the input balls.  Mostly at the precisions where the
       quadrature uses doubles; sometimes higher. */
    for (iter = 0; iter < 400 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, z, s0, z0, r, f;
        slong prec;
        int upper, ok;

        acb_init(s); acb_init(z); acb_init(s0); acb_init(z0);
        acb_init(r); acb_init(f);

        upper = n_randint(state, 2);
        prec = 10 + n_randint(state, (n_randint(state, 10) == 0) ? 500 : 190);

        gamma_rand_exact(s, state, 5);
        gamma_rand_exact(z, state, 5);
        acb_set(s0, s);
        acb_set(z0, z);

        if (iter % 2)
        {
            gamma_rand_ball(s0, s, state);
            gamma_rand_ball(z0, z, state);
        }

        gamma_reference(f, s0, z0, upper, gamma_reference_prec(s, z, prec));

        /* also with aliasing */
        if (n_randint(state, 4) == 0)
        {
            acb_set(r, z);
            ok = _acb_hypgeom_gamma_incomplete_integration(r, s, r, upper, prec);
        }
        else
            ok = _acb_hypgeom_gamma_incomplete_integration(r, s, z, upper, prec);

        if (ok && gamma_usable(f) && !acb_overlaps(r, f))
        {
            flint_printf("FAIL: enclosure\n\n");
            flint_printf("upper = %d, prec = %wd\n", upper, prec);
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("s0 = "); acb_printd(s0, 30); flint_printf("\n");
            flint_printf("z0 = "); acb_printd(z0, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
            flint_printf("f = "); acb_printd(f, 30); flint_printf("\n");
            flint_abort();
        }

        /* exact, real and positive: the contours must be found, and
           there is no cancellation */
        if (acb_is_exact(s) && acb_is_exact(z) && acb_is_real(s) && acb_is_real(z)
            && arb_is_positive(acb_realref(s)) && arb_is_positive(acb_realref(z))
            && gamma_usable(f) && (!ok || acb_rel_accuracy_bits(r) < prec - 12))
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("upper = %d, prec = %wd, ok = %d\n", upper, prec, ok);
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(s); acb_clear(z); acb_clear(s0); acb_clear(z0);
        acb_clear(r); acb_clear(f);
    }

    /* Large parameters, where the series are of no use:
       Gamma(s+1, z) = s Gamma(s, z) + z^s e^(-z) and
       gamma(s+1, z) = s gamma(s, z) - z^s e^(-z) */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, s1, z, a, b, t;
        slong prec;
        int upper, ok;

        acb_init(s); acb_init(s1); acb_init(z);
        acb_init(a); acb_init(b); acb_init(t);

        prec = 20 + n_randint(state, 150);
        upper = n_randint(state, 2);
        gamma_rand_exact(s, state, 2 + n_randint(state, 30));
        gamma_rand_exact(z, state, 2 + n_randint(state, 30));
        acb_add_ui(s1, s, 1, ARF_PREC_EXACT);

        ok = _acb_hypgeom_gamma_incomplete_integration(a, s, z, upper, prec)
          && _acb_hypgeom_gamma_incomplete_integration(b, s1, z, upper, prec);

        acb_pow(t, z, s, prec);
        acb_mul(a, a, s, prec);
        acb_neg(s1, z);
        acb_exp(s1, s1, prec);
        acb_mul(t, t, s1, prec);
        if (upper)
            acb_add(a, a, t, prec);
        else
            acb_sub(a, a, t, prec);

        if (ok && !acb_overlaps(a, b))
        {
            flint_printf("FAIL: recurrence\n\n");
            flint_printf("upper = %d, prec = %wd\n", upper, prec);
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(s); acb_clear(s1); acb_clear(z);
        acb_clear(a); acb_clear(b); acb_clear(t);
    }

    /* The bounds: |F| <= M and |F(s1, z) - F(s0, z)| <= Ds |s1 - s0| at
       two points of the input ball */
    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, z, s0, s1, f0, f1, t;
        mag_t M, Ds, lhs, rhs;
        slong prec = 100;
        int upper;

        acb_init(s); acb_init(z); acb_init(s0); acb_init(s1);
        acb_init(f0); acb_init(f1); acb_init(t);
        mag_init(M); mag_init(Ds); mag_init(lhs); mag_init(rhs);

        upper = n_randint(state, 2);
        gamma_rand_exact(s, state, 5);
        gamma_rand_exact(z, state, 5);
        acb_set(s0, s);
        gamma_rand_ball(s1, s, state);

        if (_acb_hypgeom_gamma_incomplete_integration_bounds(M, Ds, s, z, upper))
        {
            gamma_reference(f0, s0, z, upper, gamma_reference_prec(s, z, prec));
            gamma_reference(f1, s1, z, upper, gamma_reference_prec(s, z, prec));

            if (gamma_usable(f0) && gamma_usable(f1))
            {
                acb_get_mag_lower(lhs, f0);

                if (mag_cmp(lhs, M) > 0)
                {
                    flint_printf("FAIL: magnitude bound\n\n");
                    flint_printf("upper = %d\n", upper);
                    flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                    flint_printf("f0 = "); acb_printd(f0, 30); flint_printf("\n");
                    flint_printf("M = "); mag_printd(M, 10); flint_printf("\n");
                    flint_abort();
                }

                acb_sub(t, f1, f0, prec);
                acb_get_mag_lower(lhs, t);
                acb_sub(t, s1, s0, prec);
                acb_get_mag(rhs, t);
                mag_mul(rhs, rhs, Ds);

                if (mag_cmp(lhs, rhs) > 0)
                {
                    flint_printf("FAIL: derivative bound\n\n");
                    flint_printf("upper = %d\n", upper);
                    flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                    flint_printf("s1 = "); acb_printd(s1, 30); flint_printf("\n");
                    flint_printf("Ds = "); mag_printd(Ds, 10); flint_printf("\n");
                    flint_printf("lhs = "); mag_printd(lhs, 10); flint_printf("\n");
                    flint_printf("rhs = "); mag_printd(rhs, 10); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        acb_clear(s); acb_clear(z); acb_clear(s0); acb_clear(s1);
        acb_clear(f0); acb_clear(f1); acb_clear(t);
        mag_clear(M); mag_clear(Ds); mag_clear(lhs); mag_clear(rhs);
    }

    /* The automatic functions with inexact input (evaluation at the
       midpoint with derivative bounds, among other routes): compare
       with the value at a point of the input balls; often with only z
       inexact, and with balls crossing the branch cut */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, z, s0, z0, r, f;
        slong prec;
        int upper;

        acb_init(s); acb_init(z); acb_init(s0); acb_init(z0);
        acb_init(r); acb_init(f);

        upper = n_randint(state, 2);
        prec = 10 + n_randint(state, 190);
        gamma_rand_exact(s, state, 5);
        gamma_rand_exact(z, state, 5);

        if (n_randint(state, 4) == 0)
        {
            arb_zero(acb_imagref(z));
            if (n_randint(state, 2))
                arb_neg(acb_realref(z), acb_realref(z));
        }

        acb_set(s0, s);
        acb_set(z0, z);
        if (n_randint(state, 2))
            gamma_rand_ball(s0, s, state);
        gamma_rand_ball(z0, z, state);

        gamma_reference(f, s0, z0, upper, gamma_reference_prec(s, z, prec));

        if (upper)
            acb_hypgeom_gamma_upper(r, s, z, 0, prec);
        else
            acb_hypgeom_gamma_lower(r, s, z, 0, prec);

        if (gamma_usable(f) && !acb_overlaps(r, f))
        {
            flint_printf("FAIL: enclosure (automatic)\n\n");
            flint_printf("upper = %d, prec = %wd\n", upper, prec);
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("s0 = "); acb_printd(s0, 30); flint_printf("\n");
            flint_printf("z0 = "); acb_printd(z0, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
            flint_printf("f = "); acb_printd(f, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(s); acb_clear(z); acb_clear(s0); acb_clear(z0);
        acb_clear(r); acb_clear(f);
    }

    TEST_FUNCTION_END(state);
}
