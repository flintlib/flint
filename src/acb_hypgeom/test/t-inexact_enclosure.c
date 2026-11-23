/*
    Copyright (C) 2026 Fredrik Johansson

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

/* 1F1, U and 2F1 on inexact input (the value at the midpoints with a
   bound for the variation): the enclosure of a function with known
   bounds; containment of the values at points of the input boxes;
   accuracy close to that of the input in some fixed cases; the
   integrals outside of the domain of the Euler integral. */

/* f(x, y) = exp(q x) (1 + y)^5: a function with a simple bound on a
   polydisc, |f| <= exp(|q| (|Re x| + R)) (1 + |y| + R)^5 */
static void
_ie_simple_f(acb_t res, acb_srcptr w, slong prec, void * ctx)
{
    acb_t t;
    acb_init(t);
    acb_mul_si(res, w, *((slong *) ctx), prec);
    acb_exp(res, res, prec);
    acb_add_ui(t, w + 1, 1, prec);
    acb_pow_ui(t, t, 5, prec);
    acb_mul(res, res, t, prec);
    acb_clear(t);
}

static void
_ie_eval(acb_t res, int f, const acb_t a, const acb_t b, const acb_t c, const acb_t z, slong prec)
{
    if (f == 0)
        acb_hypgeom_m(res, a, b, z, 0, prec);
    else if (f == 1)
        acb_hypgeom_m(res, a, b, z, 1, prec);
    else if (f == 2)
        acb_hypgeom_u(res, a, b, z, prec);
    else if (f == 3)
        acb_hypgeom_2f1(res, a, b, c, z, 0, prec);
    else
        acb_hypgeom_2f1(res, a, b, c, z, ACB_HYPGEOM_2F1_REGULARIZED, prec);
}

/* a random point of the box x: the midpoint, a corner or inside */
static void
_ie_random_point(acb_t p, const acb_t x, flint_rand_t state)
{
    arf_t t;
    int k, mode = n_randint(state, 4);

    arf_init(t);
    acb_get_mid(p, x);

    for (k = 0; k < 2 && mode != 0; k++)
    {
        /* the radius times j / 256, j = +-256 or random */
        arf_set_mag(t, k ? arb_radref(acb_imagref(x)) : arb_radref(acb_realref(x)));
        if (mode == 3)
            arf_mul_ui(t, t, n_randint(state, 256), ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_mul_2exp_si(t, t, 8);
        if (n_randint(state, 2))
            arf_neg(t, t);
        arf_mul_2exp_si(t, t, -8);

        if (k == 0)
            arf_add(arb_midref(acb_realref(p)), arb_midref(acb_realref(p)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_add(arb_midref(acb_imagref(p)), arb_midref(acb_imagref(p)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
    }

    arf_clear(t);
}

static void
_ie_set_random(acb_t x, flint_rand_t state, double size, int cplx, int positive)
{
    double re, im = 0.0;

    re = size * (0.05 + n_randint(state, 1000) / 1000.0);
    if (!positive && n_randint(state, 3) == 0)
        re = -re;
    re += n_randint(state, 100) / 128.0;
    if (cplx)
        im = size * (n_randint(state, 2001) / 1000.0 - 1.0) * (n_randint(state, 2) ? 1.0 : 0.1);

    acb_set_d_d(x, re, im);
}

static void
_ie_add_error(acb_t x, flint_rand_t state, slong e, int cplx)
{
    mag_t m;

    mag_init(m);
    acb_get_mag(m, x);
    mag_mul_2exp_si(m, m, -e);
    if (n_randint(state, 2))
        arb_add_error_mag(acb_realref(x), m);
    if (cplx && n_randint(state, 2))
        arb_add_error_mag(acb_imagref(x), m);
    mag_clear(m);
}

/* f(z) = (1 + z)^n and its Taylor coefficients binomial(n, j)
   (1 + z)^(n - j); the crude values are too large by 2^K */
typedef struct
{
    ulong n;
    slong K;
}
_ie_binomial_ctx;

static void
_ie_binomial_coeff(acb_t res, slong j, const acb_t x, slong prec, int crude, void * vctx)
{
    _ie_binomial_ctx * ctx = (_ie_binomial_ctx *) vctx;
    arb_t c;

    if ((ulong) j > ctx->n)
    {
        acb_zero(res);
        return;
    }

    arb_init(c);
    arb_bin_uiui(c, ctx->n, j, prec + 10);
    acb_add_ui(res, x, 1, prec + 10);
    acb_pow_ui(res, res, ctx->n - j, prec + 10);
    acb_mul_arb(res, res, c, prec + 10);
    if (crude)
        acb_mul_2exp_si(res, res, ctx->K);
    arb_clear(c);
}

/* c_{j+2} = c_{j+1} (n - j - 1) / ((j + 2) (1 + m)) */
static void
_ie_binomial_recur(acb_t res, const acb_t c0, const acb_t c1, slong j, const acb_t m,
    slong prec, void * vctx)
{
    _ie_binomial_ctx * ctx = (_ie_binomial_ctx *) vctx;
    acb_t t;

    (void) c0;

    acb_init(t);
    acb_add_ui(t, m, 1, prec);
    acb_mul_si(t, t, j + 2, prec);
    acb_mul_si(res, c1, (slong) ctx->n - j - 1, prec);
    acb_div(res, res, t, prec);
    acb_clear(t);
}

/* the largest distance of the values of f at the corners of the box z
   from v (f as in _ie_eval, or (1 + z)^n for f = -1) */
static void
_ie_corner_variation(mag_t d, const acb_t v, int f, const acb_t a, const acb_t b, const acb_t c,
    const acb_t z, ulong n, slong prec)
{
    acb_t p, w;
    arf_t t;
    mag_t e;
    int k;

    acb_init(p);
    acb_init(w);
    arf_init(t);
    mag_init(e);
    mag_zero(d);

    for (k = 0; k < 4; k++)
    {
        acb_get_mid(p, z);
        arf_set_mag(t, arb_radref(acb_realref(z)));
        if (k & 1)
            arf_neg(t, t);
        arf_add(arb_midref(acb_realref(p)), arb_midref(acb_realref(p)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_mag(t, arb_radref(acb_imagref(z)));
        if (k & 2)
            arf_neg(t, t);
        arf_add(arb_midref(acb_imagref(p)), arb_midref(acb_imagref(p)), t, ARF_PREC_EXACT, ARF_RND_DOWN);

        if (f == -1)
        {
            acb_add_ui(w, p, 1, prec);
            acb_pow_ui(w, w, n, prec);
        }
        else
        {
            _ie_eval(w, f, a, b, c, p, prec);
        }

        acb_sub(w, w, v, prec);
        acb_get_mag(e, w);
        mag_max(d, d, e);
    }

    acb_clear(p);
    acb_clear(w);
    arf_clear(t);
    mag_clear(e);
}

TEST_FUNCTION_START(acb_hypgeom_inexact_enclosure, state)
{
    slong iter;

    /* the enclosure for a function with known bounds */
    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr args, pt;
        acb_t v, res, w;
        mag_t M, t;
        double R[2];
        slong q, prec, need, i;
        int status, real, crude;

        args = _acb_vec_init(2);
        pt = _acb_vec_init(2);
        acb_init(v); acb_init(res); acb_init(w);
        mag_init(M); mag_init(t);

        prec = 30 + n_randint(state, 300);
        q = 1 + n_randint(state, 50);
        real = n_randint(state, 2);
        crude = n_randint(state, 3);

        for (i = 0; i < 2; i++)
        {
            acb_randtest(args + i, state, 200, 3);
            if (real)
                arb_zero(acb_imagref(args + i));
            if (n_randint(state, 4) == 0)
                acb_get_mid(args + i, args + i);
        }

        R[0] = (1 + n_randint(state, 10)) / (double) (8 * q);
        R[1] = (1 + n_randint(state, 10)) / 16.0;

        /* M = exp(q (|Re x| + R0)) (1 + |y| + R1)^5 */
        arb_get_mag(M, acb_realref(args));
        mag_set_d(t, R[0]);
        mag_add(M, M, t);
        mag_mul_ui(M, M, q);
        mag_exp(M, M);
        acb_get_mag(t, args + 1);
        mag_add_ui(t, t, 1);
        {
            mag_t u;
            mag_init(u);
            mag_set_d(u, R[1]);
            mag_add(t, t, u);
            mag_clear(u);
        }
        mag_pow_ui(t, t, 5);
        mag_mul(M, M, t);

        /* (a crude bound: divided differences) */
        if (crude)
            mag_mul_2exp_si(M, M, n_randint(state, crude == 1 ? 20 : 200));

        acb_get_mid(pt, args);
        acb_get_mid(pt + 1, args + 1);
        _ie_simple_f(v, pt, prec, &q);

        status = _acb_hypgeom_inexact_enclosure(res, &need, v, args, 2, R, M, 0,
                    n_randint(state, 4) ? _ie_simple_f : NULL, &q, real, prec + 30);

        if (status == 0 && need != 0)
        {
            _ie_simple_f(v, pt, need, &q);
            status = _acb_hypgeom_inexact_enclosure(res, &need, v, args, 2, R, M, 0,
                        _ie_simple_f, &q, real, prec + 30);
        }

        if (status != 0)
        {
            for (i = 0; i < 4; i++)
            {
                _ie_random_point(pt, args, state);
                _ie_random_point(pt + 1, args + 1, state);
                _ie_simple_f(w, pt, prec + 100, &q);

                if (!acb_overlaps(res, w))
                {
                    flint_printf("FAIL: enclosure\n\n");
                    flint_printf("q = %wd, status = %d, R = %g, %g\n", q, status, R[0], R[1]);
                    flint_printf("x = "); acb_printd(args, 30); flint_printf("\n\n");
                    flint_printf("y = "); acb_printd(args + 1, 30); flint_printf("\n\n");
                    flint_printf("res = "); acb_printd(res, 30); flint_printf("\n\n");
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }

            /* with divided differences, close to the range: within
               2^6 of the first-order variation plus that of v */
            if (status == 2 && crude == 2 && !acb_is_exact(args) && !acb_is_exact(args + 1))
            {
                mag_t d, e;
                int small;

                mag_init(d); mag_init(e);

                /* (for boxes which are small in the polydisc, where the
                   first order is all there is: r_i < 2^-10 R_i) */
                mag_hypot(d, arb_radref(acb_realref(args)), arb_radref(acb_imagref(args)));
                mag_hypot(e, arb_radref(acb_realref(args + 1)), arb_radref(acb_imagref(args + 1)));
                small = mag_get_d(d) < ldexp(R[0], -10) && mag_get_d(e) < ldexp(R[1], -10);

                /* |f| (q r_x + 5 r_y / |1 + y|) */
                mag_hypot(d, arb_radref(acb_realref(args)), arb_radref(acb_imagref(args)));
                mag_mul_ui(d, d, q);
                mag_hypot(e, arb_radref(acb_realref(args + 1)), arb_radref(acb_imagref(args + 1)));
                mag_mul_ui(e, e, 5);
                acb_get_mid(w, args + 1);
                acb_add_ui(w, w, 1, 60);
                acb_get_mag_lower(t, w);
                mag_div(e, e, t);
                mag_add(d, d, e);
                acb_get_mag(t, v);
                mag_mul(d, d, t);
                mag_mul_2exp_si(d, d, 6);
                /* (R_1 <= 5/8 is smaller than the scale |1 + y| / 5 of
                   the second factor by more than that for large y) */
                acb_get_mag(t, w);
                if (mag_cmp_2exp_si(t, 2) > 0)
                {
                    mag_mul(d, d, t);
                    mag_mul_2exp_si(d, d, -2);
                }
                mag_hypot(e, arb_radref(acb_realref(v)), arb_radref(acb_imagref(v)));
                mag_mul_2exp_si(e, e, 6);
                mag_add(d, d, e);
                mag_hypot(e, arb_radref(acb_realref(res)), arb_radref(acb_imagref(res)));

                if (small && mag_cmp(e, d) > 0)
                {
                    flint_printf("FAIL: enclosure not close\n\n");
                    flint_printf("q = %wd, R = %g, %g\n", q, R[0], R[1]);
                    flint_printf("x = "); acb_printd(args, 30); flint_printf("\n\n");
                    flint_printf("y = "); acb_printd(args + 1, 30); flint_printf("\n\n");
                    flint_printf("v = "); acb_printd(v, 30); flint_printf("\n\n");
                    flint_printf("res = "); acb_printd(res, 30); flint_printf("\n\n");
                    flint_abort();
                }

                mag_clear(d); mag_clear(e);
            }
        }

        _acb_vec_clear(args, 2);
        _acb_vec_clear(pt, 2);
        acb_clear(v); acb_clear(res); acb_clear(w);
        mag_clear(M); mag_clear(t);
    }

    /* the bound from the absolute values of the terms */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr a, pa;
        acb_t b, z, pb, pz, w;
        mag_t M, t;
        slong p, i;
        int positive;

        a = _acb_vec_init(2);
        pa = _acb_vec_init(2);
        acb_init(b); acb_init(z); acb_init(pb); acb_init(pz); acb_init(w);
        mag_init(M); mag_init(t);

        p = 1 + n_randint(state, 2);
        positive = n_randint(state, 3) == 0;

        for (i = 0; i < 2; i++)
        {
            acb_set_d_d(a + i, (n_randint(state, 4000) - (positive ? 0.0 : 1000.0)) / 16.0 + 0.03125,
                positive ? 0.0 : (n_randint(state, 400) - 200.0) / 16.0);
            if (!positive && n_randint(state, 2))
                arb_zero(acb_imagref(a + i));
        }
        acb_set_d_d(b, (n_randint(state, 4000) - (positive ? 0.0 : 1000.0)) / 16.0 + 0.046875,
            positive ? 0.0 : (n_randint(state, 400) - 200.0) / 16.0);
        if (p == 1)
            acb_set_d_d(z, (n_randint(state, 4000) - (positive ? 0.0 : 2000.0)) / 16.0,
                positive ? 0.0 : (n_randint(state, 400) - 200.0) / 16.0);
        else
            acb_set_d_d(z, (n_randint(state, 1900) - (positive ? 0.0 : 950.0)) / 2048.0,
                positive ? 0.0 : (n_randint(state, 1200) - 600.0) / 2048.0);

        if (!positive || n_randint(state, 2))
        {
            mag_set_ui_2exp_si(t, 1 + n_randint(state, 8), -(slong) n_randint(state, 30));
            for (i = 0; i < p; i++)
                if (n_randint(state, 2))
                    acb_add_error_mag(a + i, t);
            if (n_randint(state, 2))
                acb_add_error_mag(b, t);
            mag_mul_2exp_si(t, t, -6);
            if (n_randint(state, 2))
                acb_add_error_mag(z, t);
        }

        _acb_hypgeom_inexact_terms_bound(M, a, p, b, 1, z, 20000, NULL);

        /* (with a limit for the sum: the same bound, or none if the
           sum exceeds the limit, which it does not if the limit is the
           bound itself) */
        if (mag_is_finite(M))
        {
            mag_t M2, lim;

            mag_init(M2);
            mag_init(lim);

            _acb_hypgeom_inexact_terms_bound(M2, a, p, b, 1, z, 20000, M);
            mag_mul_2exp_si(lim, M, -(slong) (iter % 4));
            _acb_hypgeom_inexact_terms_bound(lim, a, p, b, 1, z, 20000, lim);

            if (!mag_equal(M2, M) || (mag_is_finite(lim) && !mag_equal(lim, M)))
            {
                flint_printf("FAIL: bound from the terms with a limit\n\n");
                flint_printf("M = "); mag_printd(M, 10); flint_printf("\n\n");
                flint_printf("M2 = "); mag_printd(M2, 10); flint_printf("\n\n");
                flint_printf("lim = "); mag_printd(lim, 10); flint_printf("\n\n");
                flint_abort();
            }

            mag_clear(M2);
            mag_clear(lim);
        }

        if (mag_is_finite(M))
        {
            for (i = 0; i < 3; i++)
            {
                _ie_random_point(pa, a, state);
                _ie_random_point(pa + 1, a + 1, state);
                _ie_random_point(pb, b, state);
                _ie_random_point(pz, z, state);

                if (p == 1)
                    acb_hypgeom_m_1f1(w, pa, pb, pz, 0, 1000);
                else
                    acb_hypgeom_2f1_direct(w, pa, pa + 1, pb, pz, 0, 1000);

                acb_get_mag_lower(t, w);

                if (acb_is_finite(w) && mag_cmp(t, M) > 0)
                {
                    flint_printf("FAIL: bound from the terms\n\n");
                    flint_printf("p = %wd\n", p);
                    flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                    flint_printf("a2 = "); acb_printd(a + 1, 30); flint_printf("\n\n");
                    flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_printf("M = "); mag_printd(M, 10); flint_printf("\n\n");
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }

            /* positive terms and exact input: the bound is close to the
               sum (the ratios of the terms are bounded in pairs of an
               upper and a lower parameter, which costs a factor of 2
               to 3 when the parameters differ much) */
            if (positive && acb_is_exact(a) && acb_is_exact(a + 1) && acb_is_exact(b) && acb_is_exact(z))
            {
                if (p == 1)
                    acb_hypgeom_m_1f1(w, a, b, z, 0, 200);
                else
                    acb_hypgeom_2f1_direct(w, a, a + 1, b, z, 0, 200);

                acb_get_mag(t, w);
                mag_mul_2exp_si(t, t, 2);

                if (acb_is_finite(w) && mag_cmp(M, t) > 0)
                {
                    flint_printf("FAIL: bound from the terms (not close)\n\n");
                    flint_printf("p = %wd\n", p);
                    flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                    flint_printf("a2 = "); acb_printd(a + 1, 30); flint_printf("\n\n");
                    flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_printf("M = "); mag_printd(M, 10); flint_printf("\n\n");
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        _acb_vec_clear(a, 2);
        _acb_vec_clear(pa, 2);
        acb_clear(b); acb_clear(z); acb_clear(pb); acb_clear(pz); acb_clear(w);
        mag_clear(M); mag_clear(t);
    }

    /* The bound from the terms where the upper parameters must be
       paired with the lower ones in the right order, and where the
       ratios of the terms stay near |z|: 2F1(1, 1000; 2000; 0.99) =
       1.98115..., with the parameters in either order.  (With the
       pairs (1, 2000) and (1000, k + 1), the ratios are bounded by
       less than 1 only after 10^5 terms; and the tail bound
       t_0 |z| / (1 - |z|) is 99.) */
    {
        acb_ptr a;
        acb_t b, z;
        mag_t M, lo, hi;
        slong i;

        a = _acb_vec_init(2);
        acb_init(b); acb_init(z);
        mag_init(M); mag_init(lo); mag_init(hi);

        mag_set_d_lower(lo, 1.98115);
        mag_set_d(hi, 2.4);

        for (i = 0; i < 2; i++)
        {
            acb_set_ui(a + i, 1);
            acb_set_ui(a + 1 - i, 1000);
            acb_set_ui(b, 2000);
            acb_set_d(z, 0.99);

            _acb_hypgeom_inexact_terms_bound(M, a, 2, b, 1, z, 30000, NULL);

            if (mag_cmp(M, lo) < 0 || mag_cmp(M, hi) > 0)
            {
                flint_printf("FAIL: bound from the terms, 2F1(1, 1000; 2000; 0.99)\n\n");
                flint_printf("M = "); mag_printd(M, 10); flint_printf("\n\n");
                flint_abort();
            }
        }

        _acb_vec_clear(a, 2);
        acb_clear(b); acb_clear(z);
        mag_clear(M); mag_clear(lo); mag_clear(hi);
    }

    /* containment of the values at points of the input */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r, w, pa, pb, pc, pz;
        slong prec, e;
        int f, cplx, shape, i;
        double size;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r); acb_init(w);
        acb_init(pa); acb_init(pb); acb_init(pc); acb_init(pz);

        f = n_randint(state, 5);
        cplx = n_randint(state, 3) == 0;
        shape = n_randint(state, 3);
        size = (n_randint(state, 3) ? 20.0 : 150.0) * (1 + n_randint(state, 3));
        prec = 30 + n_randint(state, 150);
        if (n_randint(state, 10) == 0)
            prec += 500;
        e = 4 + n_randint(state, prec + 20);

        _ie_set_random(a, state, size, cplx, shape == 0);
        _ie_set_random(b, state, size, cplx, shape == 0);
        _ie_set_random(c, state, size, cplx, shape == 0);
        _ie_set_random(z, state, size, cplx, shape == 0 && f == 2);

        if (shape == 0)
        {
            /* the domains of the Euler integrals */
            acb_add_ui(a, a, 1, 200);
            acb_add(b, b, a, 200);
            acb_add_ui(b, b, 1, 200);
            acb_add(c, c, b, 200);
            acb_add_ui(c, c, 1, 200);
        }
        else if (shape == 1 && !cplx)
        {
            /* 1F1 and U with a < 0 < b and z > 0 */
            arb_abs(acb_realref(a), acb_realref(a));
            arb_neg(acb_realref(a), acb_realref(a));
            arb_abs(acb_realref(b), acb_realref(b));
            arb_abs(acb_realref(z), acb_realref(z));
        }

        if (f >= 3)
        {
            /* z in the unit disc or on the negative real axis */
            if (n_randint(state, 2))
                acb_set_d_d(z, n_randint(state, 1900) / 1000.0 - 0.95,
                    cplx ? n_randint(state, 1000) / 2000.0 : 0.0);
            else
                acb_set_d_d(z, -(double) n_randint(state, 30) - 0.5, cplx ? 0.25 : 0.0);
        }

        if (n_randint(state, 4))
            _ie_add_error(a, state, e, cplx);
        if (n_randint(state, 4))
            _ie_add_error(b, state, e, cplx);
        if (n_randint(state, 4))
            _ie_add_error(c, state, e, cplx);
        if (n_randint(state, 4))
            _ie_add_error(z, state, e, cplx);

        _ie_eval(r, f, a, b, c, z, prec);

        if (acb_is_finite(r))
        {
            for (i = 0; i < 3; i++)
            {
                _ie_random_point(pa, a, state);
                _ie_random_point(pb, b, state);
                _ie_random_point(pc, c, state);
                _ie_random_point(pz, z, state);

                _ie_eval(w, f, pa, pb, pc, pz, prec + 60);

                if (!acb_overlaps(r, w))
                {
                    flint_printf("FAIL: containment\n\n");
                    flint_printf("f = %d, prec = %wd\n\n", f, prec);
                    flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                    flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                    flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r); acb_clear(w);
        acb_clear(pa); acb_clear(pb); acb_clear(pc); acb_clear(pz);
    }

    /* accuracy close to that of the input (300 bits) */
    {
        static const double cases[10][5] = {
            /* f, a, b, c, z */
            { 0, 100.25, 250.5, 0, 120.75 },      /* the Euler integral */
            { 0, 100.25, 250.5, 0, -120.75 },
            { 0, -150.25, 80.5, 0, 200.75 },      /* a < 0: oscillating */
            { 1, 220.25, 90.5, 0, -150.75 },      /* (with Kummer's transformation) */
            { 2, 100.25, 250.5, 0, 120.75 },
            { 2, 100.25, -50.5, 0, 20.75 },
            { 2, -150.25, 80.5, 0, 200.75 },
            { 3, 100.25, 150.5, 320.75, 0.625 },
            { 3, 100.25, 150.5, 320.75, -7.5 },
            { 4, -60.25, 150.5, 320.75, 0.375 },
        };
        slong i;

        for (i = 0; i < 10; i++)
        {
            acb_t a, b, c, z, r;
            mag_t m;
            int k;

            acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r);
            mag_init(m);

            acb_set_d(a, cases[i][1]);
            acb_set_d(b, cases[i][2]);
            acb_set_d(c, cases[i][3]);
            acb_set_d(z, cases[i][4]);

            for (k = 0; k < 4; k++)
            {
                acb_ptr x = (k == 0) ? a : (k == 1) ? b : (k == 2) ? c : z;
                acb_get_mag(m, x);
                mag_mul_2exp_si(m, m, -300);
                arb_add_error_mag(acb_realref(x), m);
            }

            _ie_eval(r, (int) cases[i][0], a, b, c, z, 333);

            /* (the condition numbers are below 2^12) */
            if (acb_rel_accuracy_bits(r) < 300 - 12 - 12)
            {
                flint_printf("FAIL: accuracy\n\n");
                flint_printf("case %wd: accuracy %wd\n\n", i, acb_rel_accuracy_bits(r));
                flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                flint_abort();
            }

            acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r);
            mag_clear(m);
        }
    }

    /* The integrals for 1F1 and U with Re(a) < 1 (from two integrals
       for U), compared to the series */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, z, r1, r2;
        slong prec;
        int cplx, f;

        acb_init(a); acb_init(b); acb_init(z); acb_init(r1); acb_init(r2);

        prec = 30 + n_randint(state, 100);
        cplx = n_randint(state, 2);
        f = n_randint(state, 3);

        /* Re(a) < 0, Re(b - a) >= 1, Re(z) > 0 */
        acb_set_d_d(a, -(double) n_randint(state, 60) - n_randint(state, 100) / 128.0 - 0.01,
            cplx ? n_randint(state, 41) / 2.0 - 10.0 : 0.0);
        acb_set_d_d(b, n_randint(state, 80) + n_randint(state, 100) / 128.0 + 1.0,
            cplx ? n_randint(state, 41) / 2.0 - 10.0 : 0.0);
        acb_set_d_d(z, n_randint(state, 150) + n_randint(state, 100) / 128.0 + 0.5,
            cplx ? (n_randint(state, 81) - 40.0) / 2.0 : 0.0);

        if (n_randint(state, 3) == 0)
        {
            mag_t m;
            mag_init(m);
            mag_set_ui_2exp_si(m, 1, -(slong) n_randint(state, 60) - 10);
            acb_add_error_mag(a, m);
            acb_add_error_mag(z, m);
            mag_clear(m);
        }

        if (f == 2)
        {
            acb_hypgeom_u_integration(r1, a, b, z, prec);
            acb_hypgeom_u_1f1(r2, a, b, z, prec + 300);
        }
        else
        {
            acb_hypgeom_1f1_integration(r1, a, b, z, f, prec);
            acb_hypgeom_m_1f1(r2, a, b, z, f, prec + 300);
        }

        if (acb_is_finite(r1) && acb_is_finite(r2) && !acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: integration (two integrals)\n\n");
            flint_printf("f = %d, prec = %wd\n\n", f, prec);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a); acb_clear(b); acb_clear(z); acb_clear(r1); acb_clear(r2);
    }

    /* (where the function oscillates, the integrals are accurate) */
    {
        acb_t a, b, z, r1, r2;

        acb_init(a); acb_init(b); acb_init(z); acb_init(r1); acb_init(r2);

        acb_set_d(a, -244.0625);
        acb_set_d(b, 143.96875);
        acb_set_d(z, 54.25);

        acb_hypgeom_1f1_integration(r1, a, b, z, 0, 64);
        acb_hypgeom_m_1f1(r2, a, b, z, 0, 600);

        if (!acb_overlaps(r1, r2) || acb_rel_accuracy_bits(r1) < 50)
        {
            flint_printf("FAIL: integration (two integrals, accuracy)\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_hypgeom_u_integration(r1, a, b, z, 64);
        acb_hypgeom_u_1f1(r2, a, b, z, 1200);

        if (!acb_overlaps(r1, r2) || acb_rel_accuracy_bits(r1) < 50)
        {
            flint_printf("FAIL: integration (two integrals for U, accuracy)\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a); acb_clear(b); acb_clear(z); acb_clear(r1); acb_clear(r2);
    }

    /* The Taylor expansion in one inexact argument, for (1 + z)^n with
       crude bounds which are too large by 2^K: the values at points of
       the box are inside, the radius is close to the variation over
       the box, and the coefficients from the recurrence give the same
       as those from their formula. */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, r1, r2, v, w, pz;
        mag_t d, rad, tol;
        _ie_binomial_ctx ctx;
        slong prec, e, goal, order1, order2;
        int cplx, real, s1, s2, i;

        acb_init(z); acb_init(r1); acb_init(r2); acb_init(v); acb_init(w); acb_init(pz);
        mag_init(d); mag_init(rad); mag_init(tol);

        ctx.n = 1 + n_randint(state, 60);
        ctx.K = n_randint(state, 3) ? n_randint(state, 300) : 0;
        cplx = n_randint(state, 2);
        prec = 30 + n_randint(state, 200);
        e = 6 + n_randint(state, prec);

        acb_set_d_d(z, n_randint(state, 4000) / 1000.0 - 1.0,
            cplx ? n_randint(state, 4000) / 1000.0 - 2.0 : 0.0);
        _ie_add_error(z, state, e, cplx);
        real = acb_is_real(z);
        goal = FLINT_MIN(prec, acb_rel_accuracy_bits(z));
        goal = FLINT_MAX(goal, 2);

        s1 = _acb_hypgeom_inexact_taylor(r1, &order1, z, _ie_binomial_coeff, _ie_binomial_recur,
                &ctx, real, prec, goal, 80, 2);
        s2 = _acb_hypgeom_inexact_taylor(r2, &order2, z, _ie_binomial_coeff, NULL,
                &ctx, real, prec, goal, 80, 80);

        acb_get_mid(pz, z);
        acb_add_ui(v, pz, 1, prec + 100);
        acb_pow_ui(v, v, ctx.n, prec + 100);

        if (s1 != s2 || (s1 && (order1 != order2 || !acb_overlaps(r1, r2))))
        {
            flint_printf("FAIL: Taylor expansion with and without the recurrence\n\n");
            flint_printf("n = %wu, K = %wd, prec = %wd, orders %wd, %wd\n\n", ctx.n, ctx.K, prec, order1, order2);
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (s1)
        {
            for (i = 0; i < 4; i++)
            {
                _ie_random_point(pz, z, state);
                acb_add_ui(w, pz, 1, prec + 100);
                acb_pow_ui(w, w, ctx.n, prec + 100);

                if (!acb_overlaps(r1, w) || !acb_overlaps(r2, w))
                {
                    flint_printf("FAIL: Taylor expansion, containment\n\n");
                    flint_printf("n = %wu, K = %wd, prec = %wd, order %wd\n\n", ctx.n, ctx.K, prec, order1);
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }

            /* radius <= 1.3 (variation at the corners) + |f(m)| 2^(4 - goal) */
            _ie_corner_variation(d, v, -1, NULL, NULL, NULL, z, ctx.n, prec + 100);
            mag_mul_ui(d, d, 13);
            mag_div_ui(d, d, 10);
            acb_get_mag(tol, v);
            mag_mul_2exp_si(tol, tol, 4 - goal);
            mag_add(d, d, tol);
            mag_max(rad, arb_radref(acb_realref(r1)), arb_radref(acb_imagref(r1)));

            if (mag_cmp(rad, d) > 0)
            {
                flint_printf("FAIL: Taylor expansion, radius\n\n");
                flint_printf("n = %wu, K = %wd, prec = %wd, goal = %wd, order %wd\n\n", ctx.n, ctx.K, prec, goal, order1);
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("rad = "); mag_printd(rad, 10); flint_printf("\n\n");
                flint_printf("bound = "); mag_printd(d, 10); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(z); acb_clear(r1); acb_clear(r2); acb_clear(v); acb_clear(w); acb_clear(pz);
        mag_clear(d); mag_clear(rad); mag_clear(tol);
    }

    /* Polynomials (2F1 with a or c - b a nonpositive integer, 1F1 with
       a or b - a, U with a or a - b + 1) with exact parameters and an
       inexact argument: the values at points of the box are inside,
       and for input of 20 bits or more the result is within 14 bits of
       the accuracy which the input allows (estimated from the values
       at the corners) in most cases. */
    {
        slong cases = 0, good = 0;

        for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
        {
            acb_t a, b, c, z, r, v, w, pz;
            mag_t d, vm;
            slong prec, e, n, acc, att;
            int f, cplx, kind, i;

            acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r); acb_init(v); acb_init(w);
            acb_init(pz);
            mag_init(d); mag_init(vm);

            f = n_randint(state, 5);
            cplx = n_randint(state, 3) == 0;
            kind = n_randint(state, 2);
            n = 2 + n_randint(state, n_randint(state, 3) ? 60 : 300);
            prec = 30 + n_randint(state, 130);
            e = 8 + n_randint(state, prec + 10);

            /* parameters which are exact in binary */
            acb_set_d(b, (1 + n_randint(state, 400)) / 8.0);
            acb_set_d(c, (1 + n_randint(state, 400)) / 8.0);

            if (kind == 0)
            {
                /* a = -n */
                acb_set_si(a, -n);
                if (n_randint(state, 3) == 0)
                {
                    /* (Legendre) */
                    acb_set_si(b, n + 1);
                    acb_one(c);
                }
            }
            else if (f <= 1)
            {
                /* b - a = -n */
                acb_add_si(a, b, n, ARF_PREC_EXACT);
            }
            else if (f == 2)
            {
                /* a - b + 1 = -n */
                acb_set_d(a, (1 + n_randint(state, 400)) / 8.0);
                acb_add_si(b, a, n + 1, ARF_PREC_EXACT);
            }
            else
            {
                /* c - b = -n */
                acb_set_d(a, (1 + n_randint(state, 400)) / 8.0);
                acb_add_si(b, c, n, ARF_PREC_EXACT);
            }

            if (f <= 2)
                acb_set_d_d(z, (1 + n_randint(state, 2000)) / 1000.0 * n * ((f <= 1 && n_randint(state, 4) == 0) ? -1 : 1),
                    cplx ? (n_randint(state, 2000) / 1000.0 - 1.0) * n : 0.0);
            else if (kind == 0)
                acb_set_d_d(z, n_randint(state, 3000) / 1000.0 - 1.5, cplx ? n_randint(state, 2000) / 1000.0 - 1.0 : 0.0);
            else
                acb_set_d_d(z, n_randint(state, 4000) / 1000.0 - 3.5, cplx ? n_randint(state, 2000) / 1000.0 - 1.0 : 0.0);

            _ie_add_error(z, state, e, cplx);

            _ie_eval(r, f, a, b, c, z, prec);

            if (acb_is_finite(r) && !acb_is_exact(z))
            {
                for (i = 0; i < 3; i++)
                {
                    _ie_random_point(pz, z, state);
                    _ie_eval(w, f, a, b, c, pz, prec + 60);

                    if (!acb_overlaps(r, w))
                    {
                        flint_printf("FAIL: polynomial, containment\n\n");
                        flint_printf("f = %d, prec = %wd\n\n", f, prec);
                        flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                        flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                        flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                        flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                        flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                        flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                        flint_abort();
                    }
                }
            }

            /* the accuracy which the input allows, from the corners */
            if (!acb_is_exact(z) && acb_rel_accuracy_bits(z) >= 20)
            {
                acb_get_mid(pz, z);
                _ie_eval(v, f, a, b, c, pz, prec + 60);

                if (acb_is_finite(v) && acb_rel_accuracy_bits(v) >= prec)
                {
                    _ie_corner_variation(d, v, f, a, b, c, z, 0, prec + 60);
                    acb_get_mag_lower(vm, v);

                    if (mag_is_finite(d) && !mag_is_zero(d) && mag_cmp(vm, d) > 0)
                    {
                        mag_div(d, vm, d);
                        att = MAG_EXP(d);
                        att = FLINT_MIN(att, prec);

                        if (att >= 12)
                        {
                            acc = acb_is_finite(r) ? acb_rel_accuracy_bits(r) : -1;
                            cases++;
                            good += (acc >= att - 14);
                        }
                    }
                }
            }

            acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r); acb_clear(v); acb_clear(w);
            acb_clear(pz);
            mag_clear(d); mag_clear(vm);
        }

        if (cases >= 10 && good < 0.95 * cases)
        {
            flint_printf("FAIL: polynomials with an inexact argument: %wd of %wd accurate\n\n", good, cases);
            flint_abort();
        }
    }

    /* (fixed cases: the argument p/q rounded at the precision, where
       the terms cancel by 80 to 1700 bits) */
    {
        /* f, a, b, c, p, q, prec, bits */
        static const slong fixed[7][8] = {
            { 3, -50, 51, 1, 1, 3, 53, 38 },
            { 3, -1000, 1001, 1, 1, 3, 53, 36 },
            { 3, 423, 412, 153, 11, 10, 64, 44 },
            { 0, -50, 1, 0, 100, 3, 53, 45 },
            { 0, 301, 1, 0, -100, 3, 53, 40 },
            { 2, -1000, 1, 0, 31, 3, 53, 35 },
            { 2, 1000, 2000, 0, 1000, 3, 53, 35 } };
        acb_t a, b, c, z, r, w;
        slong i;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r); acb_init(w);

        for (i = 0; i < 7; i++)
        {
            acb_set_si(a, fixed[i][1]);
            acb_set_si(b, fixed[i][2]);
            acb_set_si(c, fixed[i][3]);
            acb_set_si(z, fixed[i][4]);
            acb_div_si(z, z, fixed[i][5], fixed[i][6]);

            _ie_eval(r, (int) fixed[i][0], a, b, c, z, fixed[i][6]);

            acb_get_mid(w, z);
            _ie_eval(w, (int) fixed[i][0], a, b, c, w, 300);

            if (!acb_overlaps(r, w) || acb_rel_accuracy_bits(w) < 200
                || acb_rel_accuracy_bits(r) < fixed[i][7])
            {
                flint_printf("FAIL: polynomial with a rounded argument, case %wd\n\n", i);
                flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r); acb_clear(w);
    }

    TEST_FUNCTION_END(state);
}
