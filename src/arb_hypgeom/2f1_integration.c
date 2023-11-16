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

    exp(f(t)) where f(z) = (b-1)*log(t) + (c-b-1)*log(1-t) - a*log(1-z*t)

Magnitude bound:

    |exp(f(t))| = exp(Re(f(t))) = exp(g(u,v)),  t = u+v*i

    g(u,v) = 0.5*[(b-1)*log(u^2+v^2) + (c-b-1)*log((u-1)^2+v^2) + (-a)*log((v*z)^2+(u*z-1)^2)]

Evaluating g(u,v) directly gives poor results; we get better bounds
using linearization.

    d/du g(u,v) = u*(b-1)/(u^2+v^2) + (u-1)*(c-b-1)/(v^2+(1-u)^2) + (-a)*z*(uz-1)/((vz)^2+(uz-1)^2)
    d/dv g(u,v) = v*[(b-1)/(u^2+v^2) + (c-b-1)/(v^2+(1-u)^2) + (-a)*z^2/((vz)^2+(1-uz)^2)]

*/

/* 0.5*[(b-1)*log(u^2+v^2) + (c-b-1)*log((u-1)^2+v^2) + (-a)*log((v*z)^2+(u*z-1)^2)] */
static di_t
di_integrand_edge(di_t u, di_t v, di_t b1, di_t cb1, di_t nega, di_t z)
{
    di_t X, Y, Z;

    X = di_fast_mul(b1, di_fast_log_nonnegative(di_fast_add(di_fast_sqr(u), di_fast_sqr(v))));

    if (cb1.a == 0 && cb1.b == 0)
        Y = di_interval(0.0, 0.0);
    else
        Y = di_fast_mul(cb1, di_fast_log_nonnegative(di_fast_add(di_fast_sqr(di_fast_sub_d(u, 1.0)), di_fast_sqr(v))));

    Z = di_fast_mul(nega, di_fast_log_nonnegative(di_fast_add(
        di_fast_sqr(di_fast_mul(v, z)),
        di_fast_sqr(di_fast_sub_d(di_fast_mul(u, z), 1.0)))));

    return di_fast_mul_d(di_fast_add(X, di_fast_add(Y, Z)), 0.5);
}

/*
which == 0 - d/du g(u,v) = u*(b-1)/(u^2+v^2) + (u-1)*(c-b-1)/(v^2+(1-u)^2) + (-a)*z*(uz-1)/((vz)^2+(uz-1)^2)
which == 1 - d/dv g(u,v) = v*[(b-1)/(u^2+v^2) + (c-b-1)/(v^2+(1-u)^2) + (-a)*z^2/((vz)^2+(1-uz)^2)]
*/
static di_t
di_integrand_edge_diff(di_t u, di_t v, di_t b1, di_t cb1, di_t nega, di_t z, int which)
{
    di_t X, Y, Z, uz1;

    uz1 = di_fast_sub_d(di_fast_mul(u, z), 1.0);

    X = di_fast_div(b1, di_fast_add(di_fast_sqr(u), di_fast_sqr(v)));

    if (cb1.a == 0 && cb1.b == 0)
        Y = di_interval(0.0, 0.0);
    else
        Y = di_fast_div(cb1, di_fast_add(di_fast_sqr(di_fast_sub_d(u, 1.0)), di_fast_sqr(v)));

    Z = di_fast_div(nega, di_fast_add(di_fast_sqr(di_fast_mul(v, z)), di_fast_sqr(uz1)));

    if (which == 0)
    {
        X = di_fast_mul(X, u);
        Y = di_fast_mul(Y, di_fast_sub_d(u, 1.0));
        Z = di_fast_mul(Z, di_fast_mul(z, uz1));

        return di_fast_add(X, di_fast_add(Y, Z));
    }
    else
    {
        Z = di_fast_mul(Z, di_fast_sqr(z));

        return di_fast_mul(di_fast_add(X, di_fast_add(Y, Z)), v);
    }
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
integrand_wide_bound5(acb_t res, const acb_t t, const arb_t b1, const arb_t cb1, const arb_t nega, const arb_t z, slong prec)
{
    slong i, N;
    di_t du, dv, db1, dcb1, dnega, dz, dg, dgprime;
    double radius, bound;
    double start, end;
    int which;
    arb_t abound;

    N = 8;
    bound = -D_INF;

    db1 = arb_get_di(b1);
    dcb1 = arb_get_di(cb1);
    dnega = arb_get_di(nega);
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
            dg = di_integrand_edge(du, di_fast_mid(dv), db1, dcb1, dnega, dz);
            dgprime = di_integrand_edge_diff(du, dv, db1, dcb1, dnega, dz, 1);
            dg = di_fast_add(dg, di_fast_mul(dgprime, di_interval(0.0, radius)));
#else
            dg = di_integrand_edge(du, dv, db1, dcb1, dnega, dz);
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
        dg = di_integrand_edge(di_fast_mid(du), dv, db1, dcb1, dnega, dz);
        dgprime = di_integrand_edge_diff(du, dv, db1, dcb1, dnega, dz, 0);
        dg = di_fast_add(dg, di_fast_mul(dgprime, di_interval(0.0, radius)));
#else
        dg = di_integrand_edge(du, dv, db1, dcb1, dnega, dz);
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
    arb_srcptr b1, cb1, nega, z;
    acb_t s, u, v;

    b1 = ((arb_srcptr) param) + 0;
    cb1 = ((arb_srcptr) param) + 1;
    nega = ((arb_srcptr) param) + 2;
    z = ((arb_srcptr) param) + 3;

    acb_init(s);
    acb_init(u);
    acb_init(v);

    acb_sub_ui(v, t, 1, prec);
    acb_neg(v, v);

    acb_mul_arb(u, t, z, prec);
    acb_sub_ui(u, u, 1, prec);
    acb_neg(u, u);

    if (order == 1)
    {
        if (!arb_is_positive(acb_realref(t)) || !arb_is_positive(acb_realref(u)) ||
            (!(arb_is_positive(acb_realref(v)) || arb_is_zero(cb1))))
            acb_indeterminate(out);
        else
        {
            integrand_wide_bound5(out, t, b1, cb1, nega, z, prec);

#if 0
            /* t^(b-1) */
            acb_log(s, t, prec);
            acb_mul_arb(s, s, b1, prec);
            /* (1-t)^(c-b-1) */
            acb_log(v, v, prec);
            acb_mul_arb(v, v, cb1, prec);
            /* (1-zt)^(-a) */
            acb_log(u, u, prec);
            acb_mul_arb(u, u, nega, prec);

            acb_add(out, s, u, prec);
            acb_add(out, out, v, prec);
            acb_exp(out, out, prec);
#endif
        }
    }
    else
    {
        if (acb_contains_zero(t) || acb_contains_zero(v))
        {
            /* t^(b-1) */
            acb_my_pow_arb(s, t, b1, prec);
            /* (1-t)^(c-b-1) */
            acb_my_pow_arb(v, v, cb1, prec);
            /* (1-zt)^(-a) */
            acb_my_pow_arb(u, u, nega, prec);
            acb_mul(out, s, u, prec);
            acb_mul(out, out, v, prec);
        }
        else
        {
            /* t^(b-1) */
            acb_log(s, t, prec);
            acb_mul_arb(s, s, b1, prec);
            /* (1-t)^(c-b-1) */
            acb_log(v, v, prec);
            acb_mul_arb(v, v, cb1, prec);
            /* (1-zt)^(-a) */
            if (arb_is_zero(nega))
            {
                acb_zero(u);
            }
            else
            {
                acb_log(u, u, prec);
                acb_mul_arb(u, u, nega, prec);
            }

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
estimate_magnitude(mag_t res, const arb_t ra, const arb_t rb, const arb_t rc, const arb_t rz)
{
    double a, b, c, z, t1, t2, u, m;
    fmpz_t e;

    a = arf_get_d(arb_midref(ra), ARF_RND_NEAR);
    b = arf_get_d(arb_midref(rb), ARF_RND_NEAR);
    c = arf_get_d(arb_midref(rc), ARF_RND_NEAR);
    z = arf_get_d(arb_midref(rz), ARF_RND_NEAR);

    u = 4*(b-1)*(2+a-c)*z + (2-c+(1+a-b)*z)*(2-c+(1+a-b)*z);

    if (u >= 0.0)
    {
        t1 = (2-c+z*(1+a-b) + sqrt(u)) / (2 * (2+a-c) * z);
        t2 = (2-c+z*(1+a-b) - sqrt(u)) / (2 * (2+a-c) * z);
    }
    else
    {
        t1 = 1e-8;
        t2 = 1 - 1e-8;
    }

    /* todo: more reliable solution when peak is at (or close to) 0 or 1 */
    t1 = FLINT_MAX(t1, 1e-10);
    t2 = FLINT_MAX(t2, 1e-10);
    t1 = FLINT_MIN(t1, 1 - 1e-10);
    t2 = FLINT_MIN(t2, 1 - 1e-10);

    m = -1e300;

    if (t1 > 0.0 && t1 < 1.0 && z * t1 < 1.0)
    {
        t1 = (b - 1) * log(t1) + (c - b - 1) * log(1 - t1) - a * log(1 - z * t1);
        m = FLINT_MAX(m, t1);
    }

    if (t2 > 0.0 && t2 < 1.0 && z * t2 < 1.0)
    {
        t2 = (b - 1) * log(t2) + (c - b - 1) * log(1 - t2) - a * log(1 - z * t2);
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

int
_arb_hypgeom_2f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)
{
    acb_calc_integrate_opt_t opt;
    arb_struct param[4];
    arb_t t, b1, cb1, nega;
    acb_t zero, one, I;
    mag_t abs_tol;
    int ok;

    arb_init(t);
    arb_init(b1);
    arb_init(cb1);
    arb_init(nega);

    arb_sub_ui(b1, b, 1, prec);
    arb_sub(cb1, c, b, prec);
    arb_sub_ui(cb1, cb1, 1, prec);
    arb_neg(nega, a);

    arb_one(t);
    ok = arb_is_finite(z) && arb_lt(z, t);
    ok = ok && arb_is_nonnegative(b1);
    ok = ok && arb_is_nonnegative(cb1);

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

        param[0] = *b1;
        param[1] = *cb1;
        param[2] = *nega;
        param[3] = *z;

        acb_calc_integrate_opt_init(opt);
        /* opt->verbose = 2; */
        /* opt->eval_limit = WORD_MAX; */

        acb_one(one);
        estimate_magnitude(abs_tol, a, b, c, z);
        mag_mul_2exp_si(abs_tol, abs_tol, -prec);

        acb_calc_integrate(I, integrand, param, zero, one, prec, abs_tol, opt, prec);

        if (!(regularized & 1))
        {
            arb_gamma(t, c, prec);
            arb_mul(acb_realref(I), acb_realref(I), t, prec);
        }
        arb_rgamma(t, b, prec);
        arb_mul(acb_realref(I), acb_realref(I), t, prec);
        arb_sub(t, c, b, prec);
        arb_rgamma(t, t, prec);
        arb_mul(acb_realref(I), acb_realref(I), t, prec);

        arb_set(res, acb_realref(I));

        mag_clear(abs_tol);
        acb_clear(zero);
        acb_clear(one);
        acb_clear(I);
    }

    arb_clear(t);
    arb_clear(b1);
    arb_clear(cb1);
    arb_clear(nega);

    return ok;
}

void
arb_hypgeom_2f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)
{
    arb_t res2;
    arb_init(res2);

    if (arf_cmp(arb_midref(a), arb_midref(b)) < 0)
    {
        if (!_arb_hypgeom_2f1_integration(res2, a, b, c, z, regularized, prec))
            _arb_hypgeom_2f1_integration(res2, b, a, c, z, regularized, prec);
    }
    else
    {
        if (!_arb_hypgeom_2f1_integration(res2, b, a, c, z, regularized, prec))
            _arb_hypgeom_2f1_integration(res2, a, b, c, z, regularized, prec);
    }

    arb_swap(res, res2);
    arb_clear(res2);
}
