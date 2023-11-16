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

Integrand (see comments for 1f1_integration):

    exp(f(t)) where f(z) = -z*t + a1*log(t) + ba1*log(1+t)

    g(u,v) = -z*u + 0.5*[a1*log(u^2+v^2) + ba1*log((1+u)^2+v^2)]

    d/du g(u,v) = -z + u*a1/(u^2+v^2) + (1+u)*ba1/(v^2+(1+u)^2)
    d/dv g(u,v) = v*a1/(u^2+v^2) + v*ba1/(v^2+(1+u)^2)
*/

static di_t
di_integrand_edge(di_t u, di_t v, di_t a1, di_t ba1, di_t z)
{
    di_t X, Y, Z;

    X = di_neg(di_fast_mul(z, u));

    if (a1.a == 0.0 && a1.b == 0.0)
        Y = di_interval(0.0, 0.0);
    else
        Y = di_fast_mul(a1, di_fast_log_nonnegative(di_fast_add(di_fast_sqr(u), di_fast_sqr(v))));

    Z = di_fast_mul(ba1, di_fast_log_nonnegative(di_fast_add(di_fast_sqr(di_fast_add_d(u, 1.0)), di_fast_sqr(v))));

    return di_fast_add(X, di_fast_mul_d(di_fast_add(Y, Z), 0.5));
}

/*
which == 0 - d/du g(u,v) = -z + u*(a-1)/(u^2+v^2) + (u+1)*(b-a-1)/(v^2+(1+u)^2)
which == 1 - d/dv g(u,v) = v*(a-1)/(u^2+v^2) + v*(b-a-1)/(v^2+(1+u)^2)
*/
static di_t
di_integrand_edge_diff(di_t u, di_t v, di_t a1, di_t ba1, di_t z, int which)
{
    di_t Y, Z;

    if (a1.a == 0.0 && a1.b == 0.0)
        Y = di_interval(0.0, 0.0);
    else
        Y = di_fast_div(a1, di_fast_add(di_fast_sqr(u), di_fast_sqr(v)));

    Z = di_fast_div(ba1, di_fast_add(di_fast_sqr(di_fast_add_d(u, 1.0)), di_fast_sqr(v)));

    if (which == 0)
        return di_fast_add(di_neg(z), di_fast_add(di_fast_mul(u, Y), di_fast_mul(di_fast_add_d(u, 1.0), Z)));
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

    acb_add_ui(v, t, 1, prec);

    if (order == 1)
    {
        if (!(arb_is_positive(acb_realref(t)) || arb_is_zero(a1)) ||
            !arb_is_positive(acb_realref(v)))
            acb_indeterminate(out);
        else
            integrand_wide_bound5(out, t, a1, ba1, z, prec);
    }
    else
    {
        if (acb_contains_zero(t) || acb_contains_zero(v))
        {
            /* exp(-z t) */
            acb_mul_arb(s, t, z, prec);
            acb_neg(s, s);
            acb_exp(s, s, prec);

            /* t^(a-1) */
            acb_my_pow_arb(u, t, a1, prec);

            /* (1+t)^(b-a-1) */
            acb_pow_arb(v, v, ba1, prec);

            acb_mul(out, s, u, prec);
            acb_mul(out, out, v, prec);
        }
        else
        {
            acb_mul_arb(s, t, z, prec);
            acb_neg(s, s);

            /* t^(a-1) */
            if (arb_is_zero(a1))
            {
                acb_zero(u);
            }
            else
            {
                acb_log(u, t, prec);
                acb_mul_arb(u, u, a1, prec);
            }

            /* (1+t)^(b-a-1) */
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
        t1 = (-2 + b - z + sqrt(u)) / (2 * z);
        t2 = (-2 + b - z - sqrt(u)) / (2 * z);
    }
    else
    {
        t1 = 1e-8;
        t2 = 1e-8;
    }

    /* todo: better estimate when peak is at 0 */
    t1 = FLINT_MAX(t1, 1e-8);
    t2 = FLINT_MAX(t2, 1e-8);

    m = -1e300;

    if (t1 > 0.0)
    {
        t1 = -z * t1 + (a - 1) * log(t1) + (b - a - 1) * log(1 + t1);
        m = FLINT_MAX(m, t1);
    }

    if (t2 > 0.0)
    {
        t2 = -z * t2 + (a - 1) * log(t2) + (b - a - 1) * log(1 + t2);
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

static void
bound_tail(mag_t bound, const arb_t a1, const arb_t ba1, const arb_t z, const arb_t N, slong prec)
{
    arb_t s, u, v, C;

    arb_init(s);
    arb_init(u);
    arb_init(v);
    arb_init(C);

    /*
    Assume N >= 1 and t >= 0.

        -z*(N+t) + (a-1)*log(N+t) + (b-a-1)*log(1+N+t)
    <=  [-z*N + (a-1)*log(N) + (b-a-1)*log(1+N)] + [-z*t + [max(0, a-1) + max(0, b-a-1)]*log(1+t/N)]
    <=  [-z*N + (a-1)*log(N) + (b-a-1)*log(1+N)] + [-z*t + [max(0, a-1) + max(0, b-a-1)]*(t/N)]

    Let C = max(0, a-1) + max(0, b-a-1). Then the remainder integral
    is bounded by integrand(N) * N / (N*z - C), assuming that N*z > C.
    */

    arb_max(u, u, a1, prec);
    arb_max(v, v, ba1, prec);
    arb_add(C, u, v, prec);
    /* s = N*z - C */
    arb_mul(s, N, z, prec);
    arb_sub(s, s, C, prec);

    if (arb_is_positive(s))
    {
        arb_div(C, N, s, prec);

        /* exp(-z*N) */
        arb_mul(s, N, z, prec);
        arb_neg(s, s);

        /* N^(a-1) */
        arb_log(u, N, prec);
        arb_mul(u, u, a1, prec);

        /* (1+N)^(b-a-1) */
        arb_add_ui(v, N, 1, prec);
        arb_log(v, v, prec);
        arb_mul(v, v, ba1, prec);

        arb_add(s, s, u, prec);
        arb_add(s, s, v, prec);
        arb_exp(s, s, prec);

        arb_mul(s, s, C, prec);

        arb_get_mag(bound, s);
    }
    else
    {
        mag_inf(bound);
    }

    arb_clear(s);
    arb_clear(u);
    arb_clear(v);
    arb_clear(C);
}

int
_arb_hypgeom_u_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    acb_calc_integrate_opt_t opt;
    arb_struct param[3];
    arb_t t, a1, ba1;
    acb_t zero, N, I;
    mag_t abs_tol, tail_bound;
    slong i;
    fmpz_t n;
    int ok;

    arb_init(t);
    arb_init(a1);
    arb_init(ba1);

    arb_sub_ui(a1, a, 1, prec);
    arb_sub(ba1, b, a, prec);
    arb_sub_ui(ba1, ba1, 1, prec);

    ok = arb_is_finite(z) && arb_is_positive(z);
    ok = ok && arb_is_nonnegative(a1);
    ok = ok && arb_is_finite(b);

    if (!ok)
    {
        arb_indeterminate(res);
    }
    else
    {
        mag_init(abs_tol);
        mag_init(tail_bound);
        acb_init(zero);
        acb_init(zero);
        acb_init(N);
        acb_init(I);
        fmpz_init(n);

        param[0] = *a1;
        param[1] = *ba1;
        param[2] = *z;

        acb_calc_integrate_opt_init(opt);
        /* opt->verbose = 2; */
        /* opt->eval_limit = WORD_MAX; */

        estimate_magnitude(abs_tol, a, b, z);
        mag_mul_2exp_si(abs_tol, abs_tol, -prec);

        for (i = 1; i < FLINT_BITS; i++)
        {
            fmpz_one(n);
            fmpz_mul_2exp(n, n, i);
            acb_one(N);
            arb_mul_2exp_fmpz(acb_realref(N), acb_realref(N), n);
            bound_tail(tail_bound, a1, ba1, z, acb_realref(N), 64);
            if (mag_cmp(tail_bound, abs_tol) < 0)
                break;
        }

        acb_calc_integrate(I, integrand, param, zero, N, prec, abs_tol, opt, prec);
        arb_add_error_mag(acb_realref(I), tail_bound);

        arb_rgamma(t, a, prec);
        arb_mul(acb_realref(I), acb_realref(I), t, prec);

        arb_set(res, acb_realref(I));

        mag_clear(abs_tol);
        mag_clear(tail_bound);
        acb_clear(zero);
        acb_clear(N);
        acb_clear(I);
        fmpz_clear(n);
    }

    arb_clear(t);
    arb_clear(a1);
    arb_clear(ba1);

    return ok;
}

void
arb_hypgeom_u_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    arb_t res2;
    arb_init(res2);

    if (!_arb_hypgeom_u_integration(res2, a, b, z, prec))
    {
        arb_t c, d;
        arb_init(c);
        arb_init(d);
        arb_sub(c, a, b, prec);
        arb_add_ui(c, c, 1, prec);
        arb_sub_ui(d, b, 2, prec);
        arb_neg(d, d);

        if (_arb_hypgeom_u_integration(res2, c, d, z, prec))
        {
            arb_sub_ui(c, b, 1, prec);
            arb_neg(c, c);
            arb_pow(c, z, c, prec);
            arb_mul(res2, res2, c, prec);
        }

        arb_clear(c);
        arb_clear(d);
    }

    arb_swap(res, res2);
    arb_clear(res2);
}
