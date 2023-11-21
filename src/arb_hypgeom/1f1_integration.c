/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_calc.h"
#include "double_interval.h"

/*

Integrand:

    exp(f(t)) where f(z) = z*t + (a-1)*log(t) + (b-a-1)*log(1-t)

Magnitude bound:

    |exp(f(t))| = exp(Re(f(t))) = exp(g(u,v)),  t = u+v*i

    g(u,v) = z*u + 0.5*[(a-1)*log(u^2+v^2) + (b-a-1)*log((u-1)^2+v^2)]

Evaluating g(u,v) directly gives poor results; we get better bounds
using linearization.

    d/du g(u,v) = z + u*(a-1)/(u^2+v^2) + (u-1)*(b-a-1)/(v^2+(1-u)^2)
    d/dv g(u,v) = v*(a-1)/(u^2+v^2) + v*(b-a-1)/(v^2+(1-u)^2)

Finding the extrema of g(u,v) is doable (solutions of a degree-4
polynomial) but for simplicity we just do interval arithmetic.

*/

/* z*u + 0.5 [(a-1) log(u^2+v^2) + (b-a-1) log((u-1)^2+v^2)] */
static di_t
di_integrand_edge(di_t u, di_t v, di_t a1, di_t ba1, di_t z)
{
    di_t X, Y, Z;

    X = di_fast_mul(z, u);
    Y = di_fast_mul(a1, di_fast_log_nonnegative(di_fast_add(di_fast_sqr(u), di_fast_sqr(v))));
    Z = di_fast_mul(ba1, di_fast_log_nonnegative(di_fast_add(di_fast_sqr(di_fast_sub_d(u, 1.0)), di_fast_sqr(v))));

    return di_fast_add(X, di_fast_mul_d(di_fast_add(Y, Z), 0.5));
}

/*
which == 0 - d/du g(u,v) = z + u*(a-1)/(u^2+v^2) + (u-1)*(b-a-1)/(v^2+(1-u)^2)
which == 1 - d/dv g(u,v) = v*(a-1)/(u^2+v^2) + v*(b-a-1)/(v^2+(1-u)^2)
*/
static di_t
di_integrand_edge_diff(di_t u, di_t v, di_t a1, di_t ba1, di_t z, int which)
{
    di_t Y, Z;

    Y = di_fast_div(a1, di_fast_add(di_fast_sqr(u), di_fast_sqr(v)));
    Z = di_fast_div(ba1, di_fast_add(di_fast_sqr(di_fast_sub_d(u, 1.0)), di_fast_sqr(v)));

    if (which == 0)
        return di_fast_add(z, di_fast_add(di_fast_mul(u, Y), di_fast_mul(di_fast_sub_d(u, 1.0), Z)));
    else
        return di_fast_mul(v, di_fast_add(Y, Z));
}

static di_t di_subinterval(di_t x, slong i, slong N)
{
    di_t res;
    double step;

    step = (x.b - x.a) / N;

    res.a = x.a + step * i;
    res.b = (i == N - 1) ? x.b : x.a + step * (i + 1);

    return res;
}

static void
integrand_wide_bound5(acb_t res, const acb_t t, const arb_t a1, const arb_t ba1, const arb_t z, slong prec)
{
    slong i, N;
    di_t du, dv, da1, dba1, dz, dg, dgprime;
    double radius, bound;
    double start, end;
    int which;
    arb_t abound;

    N = 8;
    bound = -D_INF;

    da1 = arb_get_di(a1);
    dba1 = arb_get_di(ba1);
    dz = arb_get_di(z);

    /* left edge: left(u) + [0, right(v)] */
    /* right edge: right(u) + [0, right(v)] */
    for (which = 0; which < 2; which++)
    {
        du = arb_get_di(acb_realref(t));
        if (which == 0)
            du.b = du.a;
        else
            du.a = du.b;

        dv = arb_get_di(acb_imagref(t));
        start = 0.0;
        end = dv.b;

        for (i = 0; i < N; i++)
        {
            dv = di_subinterval(di_interval(start, end), i, N);
            radius = di_fast_ubound_radius(dv);

            /* g(u,mid(v)) + g'(u,v) * [0, radius] */
#if 1
            dg = di_integrand_edge(du, di_fast_mid(dv), da1, dba1, dz);
            dgprime = di_integrand_edge_diff(du, dv, da1, dba1, dz, 1);
            dg = di_fast_add(dg, di_fast_mul(dgprime, di_interval(0.0, radius)));
#else
            dg = di_integrand_edge(du, dv, da1, dba1, dz);
#endif

            bound = FLINT_MAX(bound, dg.b);
        }
    }

    du = arb_get_di(acb_realref(t));
    start = du.a;
    end = du.b;

    dv = arb_get_di(acb_imagref(t));
    dv.a = dv.b;

    /* top edge: [left(u), right(u)] + right(v) */
    for (i = 0; i < N; i++)
    {
        du = di_subinterval(di_interval(start, end), i, N);
        radius = di_fast_ubound_radius(du);

        /* g(mid(u),v) + g'(u,v) * [0, radius] */
#if 1
        dg = di_integrand_edge(di_fast_mid(du), dv, da1, dba1, dz);
        dgprime = di_integrand_edge_diff(du, dv, da1, dba1, dz, 0);
        dg = di_fast_add(dg, di_fast_mul(dgprime, di_interval(0.0, radius)));
#else
        dg = di_integrand_edge(du, dv, da1, dba1, dz);
#endif

        bound = FLINT_MAX(bound, dg.b);
    }

    arb_init(abound);
    arb_set_d(abound, bound);
    arb_exp(abound, abound, prec);

    acb_zero(res);
    arb_add_error(acb_realref(res), abound);
    arb_add_error(acb_imagref(res), abound);

    arb_clear(abound);
}

/* todo: fix acb_pow(_arb) */
static void
acb_my_pow_arb(acb_t res, const acb_t a, const arb_t b, slong prec)
{
    if (acb_contains_zero(a) && arb_is_positive(b))
    {
        /* |a^b| <= |a|^b */
        arb_t t, u;

        arb_init(t);
        arb_init(u);

        acb_abs(t, a, prec);
        arb_get_abs_ubound_arf(arb_midref(t), t, prec);
        mag_zero(arb_radref(t));

        if (arf_cmpabs_2exp_si(arb_midref(t), 0) < 0)
            arb_get_abs_lbound_arf(arb_midref(u), b, prec);
        else
            arb_get_abs_ubound_arf(arb_midref(u), b, prec);

        arb_pow(t, t, u, prec);

        acb_zero(res);
        acb_add_error_arb(res, t);

        arb_clear(t);
        arb_clear(u);
    }
    else
    {
        acb_pow_arb(res, a, b, prec);
    }
}

static int
integrand(acb_ptr out, const acb_t t, void * param, slong order, slong prec)
{
    arb_srcptr a1, ba1, z;
    acb_t s, u, v;

    a1 = ((arb_srcptr) param) + 0;
    ba1 = ((arb_srcptr) param) + 1;
    z = ((arb_srcptr) param) + 2;

    acb_init(s);
    acb_init(u);
    acb_init(v);

    acb_sub_ui(v, t, 1, prec);
    acb_neg(v, v);

    if (order == 1)
    {
        if (!arb_is_positive(acb_realref(t)) || !arb_is_positive(acb_realref(v)))
            acb_indeterminate(out);
        else
            integrand_wide_bound5(out, t, a1, ba1, z, prec);
    }
    else
    {
        if (acb_contains_zero(t) || acb_contains_zero(v))
        {
            /* exp(z t) */
            acb_mul_arb(s, t, z, prec);
            acb_exp(s, s, prec);

            /* t^(a-1) */
            acb_my_pow_arb(u, t, a1, prec);

            /* (1-t)^(b-a-1) */
            acb_my_pow_arb(v, v, ba1, prec);

            acb_mul(out, s, u, prec);
            acb_mul(out, out, v, prec);
        }
        else
        {
            acb_mul_arb(s, t, z, prec);

            /* t^(a-1) */
            acb_log(u, t, prec);
            acb_mul_arb(u, u, a1, prec);

            /* (1-t)^(b-a-1) */
            acb_log(v, v, prec);
            acb_mul_arb(v, v, ba1, prec);

            acb_add(out, s, u, prec);
            acb_add(out, out, v, prec);
            acb_exp(out, out, prec);
        }
    }

    acb_clear(s);
    acb_clear(u);
    acb_clear(v);

    return 0;
}

/* estimate integral by magnitude at peak */
static void
estimate_magnitude(mag_t res, const arb_t ra, const arb_t rb, const arb_t rz)
{
    double a, b, z, t1, t2, u, m;
    fmpz_t e;

    a = arf_get_d(arb_midref(ra), ARF_RND_NEAR);
    b = arf_get_d(arb_midref(rb), ARF_RND_NEAR);
    z = arf_get_d(arb_midref(rz), ARF_RND_NEAR);

    u = 4 - 4*b + b*b + 4*a*z - 2*b*z + z*z;

    if (u >= 0.0)
    {
        t1 = (2 - b + z + sqrt(u)) / (2 * z);
        t2 = (2 - b + z - sqrt(u)) / (2 * z);
    }
    else
    {
        t1 = 1e-8;
        t2 = 1 - 1e-8;
    }

    m = -1e300;

    if (t1 > 0.0 && t1 < 1.0)
    {
        t1 = z * t1 + (a - 1) * log(t1) + (b - a - 1) * log(1 - t1);
        m = FLINT_MAX(m, t1);
    }

    if (t2 > 0.0 && t2 < 1.0)
    {
        t2 = z * t2 + (a - 1) * log(t2) + (b - a - 1) * log(1 - t2);
        m = FLINT_MAX(m, t2);
    }

    m /= log(2);

    if (fabs(m) < 1e300)
    {
        fmpz_init(e);
        fmpz_set_d(e, m);
        mag_set_d_2exp_fmpz(res, 1.0, e);
        fmpz_clear(e);
    }
    else
    {
        mag_zero(res);
    }
}

void
arb_hypgeom_1f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)
{
    acb_calc_integrate_opt_t opt;
    arb_struct param[3];
    arb_t t, a1, ba1;
    acb_t zero, one, I;
    mag_t abs_tol;
    int ok;

    arb_init(t);
    arb_init(a1);
    arb_init(ba1);

    arb_sub_ui(a1, a, 1, prec);
    arb_sub(ba1, b, a, prec);
    arb_sub_ui(ba1, ba1, 1, prec);

    ok = arb_is_finite(z);
    ok = ok && arb_is_nonnegative(a1);
    ok = ok && arb_is_nonnegative(ba1);

    if (!ok)
    {
        arb_indeterminate(res);
    }
    else
    {
        mag_init(abs_tol);
        acb_init(zero);
        acb_init(one);
        acb_init(I);

        param[0] = *a1;
        param[1] = *ba1;
        param[2] = *z;

        acb_calc_integrate_opt_init(opt);
        /* opt->verbose = 2; */
        /* opt->eval_limit = WORD_MAX; */

        acb_one(one);
        estimate_magnitude(abs_tol, a, b, z);
        mag_mul_2exp_si(abs_tol, abs_tol, -prec);
        acb_calc_integrate(I, integrand, param, zero, one, prec, abs_tol, opt, prec);

        if (!regularized)
        {
            arb_gamma(t, b, prec);
            arb_mul(acb_realref(I), acb_realref(I), t, prec);
        }
        arb_rgamma(t, a, prec);
        arb_mul(acb_realref(I), acb_realref(I), t, prec);
        arb_sub(t, b, a, prec);
        arb_rgamma(t, t, prec);
        arb_mul(acb_realref(I), acb_realref(I), t, prec);

        arb_set(res, acb_realref(I));

        mag_clear(abs_tol);
        acb_clear(zero);
        acb_clear(one);
        acb_clear(I);
    }

    arb_clear(t);
    arb_clear(a1);
    arb_clear(ba1);
}
