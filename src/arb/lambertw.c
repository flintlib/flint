/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "double_extras.h"

/* Helper functions to compute W_{-1}(x) on (-1/e,0) in double precision --
   just to get a good starting value for the multiprecision code, and
   not optimized for accuracy. Has underflow problems very close to 0.
   In part based on d_lambertw in flint/double_extras which implements
   the principal branch. */

#define ONE_OVER_E ldexp(6627126856707895.0, -54)
#define CORRECTION 4.3082397558469466e-17
#define POLY(p, x) d_polyval((p), sizeof(p) / sizeof(double), (x))

static double
d_halley(double x, double w)
{
    double t, u, v;
    t = exp(w);
    u = 2*w + 2;
    v = w*t - x;
    t = w - u*v / (u*t*(w+1) - (w+2)*v);
    return t;
}

static const double pol1[11] = {
    -1.0000000000000000000, 2.3316439815971242034, -1.8121878856393634902,
    1.9366311144923597554, -2.3535512018816145168, 3.0668589010506319129,
    -4.1753356002581771389, 5.8580237298747741488, -8.4010322175239773710,
    12.250753501314460424, -18.100697012472442755 };

static const double pol2[4] = {
    -5.2012020327515463962,-24.075862656446909233,
    -26.500221957196285366,2.3340178581744999812 };

static const double pol3[4] = {
    1.0000000000000000000,0.14831080741950550754,
    -13.649088405005569258,-18.975103873227202378 };

static const double pol4[5] = {
    -8.4834127832006526854,634.84191267691313548,-2640.6635889188399862,
    -12935.640726994524734,-7875.3418281832619890 };

static const double pol5[5] = {
    1.0,-121.07185283214167442,1287.5430771188798866,
    1550.0693150055578327,-3278.4808321541988073 };

static const double pol6[5] = {
    -12.169991898228748602,32778.761570863291802,-1.0480461503378867869e7,
    4.7898751364090879209e8,-7.802332913704000874e8 };

static const double pol7[5] = {
    1.0,-3556.4306263369027831,1.4761527435056145298e6,
    -9.8425904825010893103e7,7.0373606710750560344e8 };

double
d_lambertw_branch1(double x)
{
    double w, u;

    if (x < -ONE_OVER_E || x >= 0.0)
    {
        return D_NAN;
    }
    else if (x < -ONE_OVER_E + 1/32.)
    {
        w = POLY(pol1, -sqrt((x + ONE_OVER_E) + CORRECTION));
        if (x + ONE_OVER_E > 0.0003)
            w = d_halley(x, w);
        return w;
    }
    else if (x <= -1/4.)
    {
        return d_halley(x, POLY(pol2,x) / POLY(pol3,x));
    }
    else if (x < -1/32.)
    {
        return d_halley(x, POLY(pol4,x) / POLY(pol5,x));
    }
    else if (x < - 1/1024.)
    {
        return d_halley(x, POLY(pol6,x) / POLY(pol7,x));
    }
    else
    {
        w = log(-x);
        u = log(-w);
        w = (2*w*w*w - 2*(1+(w-1)*w)*u + u*u)/(2*w*w);
        if (x < -1e-15)
            w = d_halley(x, w);
        return d_halley(x, w);
    }
}

/* If branch == 0, bounds |W'(x)|. If branch != 0, bounds W_{-1}'(x).
   For the principal branch:
       For x >= 0,    W'(x) <= 1/(x+1).
       For x > -1/e,  W'(x) <  2 / sqrt(1+e*x).
    For the -1 branch:
        |W'(x)| <= 2 / sqrt(1+e*x) + 2/|x|.
*/
void
arb_lambertw_bound_prime(mag_t w, const arb_t x, int branch, slong prec)
{
    arb_t t;
    mag_t u, v;

    arb_init(t);
    mag_init(u);
    mag_init(v);

    if (arb_is_nonnegative(x) && branch == 0)
    {
        arb_get_mag_lower(w, x);
        mag_one(u);
        mag_add_lower(w, w, u);
        mag_div(w, u, w);
    }
    else
    {
        arb_const_e(t, prec);
        arb_mul(t, t, x, prec);
        arb_add_ui(t, t, 1, prec);
        arb_get_mag_lower(w, t);
        mag_rsqrt(w, w);
        mag_mul_2exp_si(w, w, 1);

        if (branch != 0)
        {
            if (arb_is_negative(x))
            {
                arb_get_mag_lower(u, x);
                mag_set_ui(v, 2);
                mag_div(v, v, u);
                mag_add(w, w, v);
            }
            else
            {
                mag_inf(w);
            }
        }
    }

    arb_clear(t);
    mag_clear(u);
    mag_clear(v);
}

/* Given an approximation w for W(x), compute a rigorous error bound.
   The precomputed value ew = e^w is optional. */
void
arb_lambertw_bound_error(mag_t res, const arb_t x, const arf_t w,
        const arb_t ew, int branch, slong prec)
{
    arb_t r, x2;
    mag_t m;

    /* Make sure that we are somewhere on the right branch. */
    if ((branch == 0 && arf_cmp_si(w, -1) < 0) ||
        (branch == 1 && arf_cmp_si(w, -1) > 0))
    {
        mag_inf(res);
        return;
    }

    arb_init(r);
    arb_init(x2);
    mag_init(m);

    if (ew != NULL)
    {
        arb_set(r, ew);
    }
    else
    {
        arb_set_arf(r, w);
        arb_exp(r, r, prec);
    }

    /* x2 = w e^w */
    arb_mul_arf(x2, r, w, prec);
    /* r = x2 - x */
    arb_sub(r, x2, x, prec);

    arb_get_mag(m, r);

    /* x2 = min(x, x+r) (W'(t) is decreasing with t) */
    if (branch == 0)
    {
        arb_min(x2, x, x2, prec);
    }
    else
    {
        arb_union(x2, x, x2, prec);
    }

    arb_lambertw_bound_prime(res, x2, branch, prec);
    mag_mul(res, res, m);

    arb_clear(r);
    arb_clear(x2);
    mag_clear(m);
}

/*
Halley iteration:
    res =  w - (w e^w - x) / (e^w (w+1) - (w+2)(w e^w - x) / (2w+2))

If certify is set, writes rigorous error bound to the radius.
*/
static void
arb_lambertw_halley_step(arb_t res, const arb_t x, const arf_t w,
                                    int branch, int certify, slong prec)
{
    arf_t t, u, v;
    arb_t ew;
    mag_t err;

    arb_init(ew);
    arf_init(t);
    arf_init(u);
    arf_init(v);
    mag_init(err);

    arb_set_arf(ew, w);
    arb_exp(ew, ew, prec);   /* todo: extra precision with large exponents? */

    arf_add_ui(u, w, 2, prec, ARF_RND_DOWN);
    arf_add_ui(v, w, 1, prec, ARF_RND_DOWN);
    arf_mul_2exp_si(v, v, 1);
    arf_div(v, u, v, prec, ARF_RND_DOWN);                 /* v = (w + 2) / (2w + 2) */

    arf_mul(t, arb_midref(ew), w, prec, ARF_RND_DOWN);    /* t = w e^w */
    arf_sub(u, t, arb_midref(x), prec, ARF_RND_DOWN);     /* u = w e^w - x */

    arf_mul(v, v, u, prec, ARF_RND_DOWN);
    arf_neg(v, v);
    arf_add(v, v, t, prec, ARF_RND_DOWN);
    arf_add(v, v, arb_midref(ew), prec, ARF_RND_DOWN);

    arf_div(t, u, v, prec, ARF_RND_DOWN);
    /* t is our new approximation for W */
    arf_sub(t, w, t, prec, ARF_RND_DOWN);

    if (certify)
    {
        arb_t et;
        arb_init(et);
        /* Inverse: x2 = t e^t. We already have e^w, so
           compute e^t = e^w e^(t-w). */
        arb_set_arf(et, w);
        arb_sub_arf(et, et, t, prec);
        arb_neg(et, et);
        arb_exp(et, et, prec);
        arb_mul(et, et, ew, prec);
        arb_lambertw_bound_error(err, x, t, et, branch, prec);
        arb_clear(et);
    }

    arf_swap(arb_midref(res), t);
    mag_swap(arb_radref(res), err);

    arb_clear(ew);
    arf_clear(t);
    arf_clear(u);
    arf_clear(v);
    mag_clear(err);
}

/* Double precision approximation good for x >= 2^1000, or
   roughly |x| <= 2^(1000-60) for the -1 branch. */
slong
arb_lambertw_initial_asymp1(arf_t res, const arf_t x, int branch, slong prec)
{
    fmpz_t e;
    double l, ll, h, t2, t3, t4;

    fmpz_init(e);
    arf_frexp(res, e, x);

    l = arf_get_d(res, ARF_RND_DOWN);
    if (branch)
        l = -l;
    l = log(l);
    l += fmpz_get_d(e) * 0.6931471805599453;

    ll = l;
    if (branch)
        ll = -ll;
    ll = log(ll);

    h = 1.0 / l;
    t2 = ll * (ll - 2) * 0.5;
    t3 = ll * (6 + ll * (2 * ll - 9)) * (1.0 / 6.0);
    t4 = ll * (-12 + ll * (36 + ll*(-22 + 3*ll))) * (1.0 / 12.0);
    l = l - ll + h * (ll + h * (t2 + h * (t3 + t4 * h)));

    arf_set_d(res, l);
    fmpz_clear(e);

    return 50;
}

/* First terms of asymptotic series, with higher precision, for huge
   exponents... */
static void
_arf_log(arf_t res, const arf_t x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_arf(t, x);
    arb_log(t, t, prec);
    arf_swap(res, arb_midref(t));
    arb_clear(t);
}

slong
arb_lambertw_initial_asymp2(arf_t res, const arf_t x, int branch, slong prec)
{
    arf_t l, ll;
    slong wp, acc;

    acc = 2 * fmpz_bits(ARF_EXPREF(x)) - 10;

    arf_init(l);
    arf_init(ll);

    wp = acc + 4;

    if (branch)
    {
        arf_neg(l, x);
        _arf_log(l, l, wp);
        arf_neg(ll, l);
        _arf_log(ll, ll, wp);
    }
    else
    {
        _arf_log(l, x, wp);
        _arf_log(ll, l, wp);
    }

    arf_div(res, ll, l, wp, ARF_RND_DOWN);
    arf_sub(res, res, ll, wp, ARF_RND_DOWN);
    arf_add(res, res, l, wp, ARF_RND_DOWN);

    arf_clear(l);
    arf_clear(ll);

    return acc;
}

/*
Computes initial approximation of W(x). Returns estimated accuracy
in bits, clamped between 0 and a reasonable value related to
the bit length of the exponent of x; thus the return value plus 2
can be used as a precision.
*/
slong
arb_lambertw_initial(arf_t res, const arf_t x, int branch, slong prec)
{
    if (arf_cmp_d(x, -ONE_OVER_E + 0.001) >= 0)
    {
        if (branch == 0)
        {
            if (arf_cmpabs_2exp_si(x, -prec) < 0)
            {
                arf_set(res, x);
                return prec;
            }
            else if (arf_cmpabs_2exp_si(x, -30) < 0)
            {
                slong acc;
                arf_set(res, x);
                acc = -arf_abs_bound_lt_2exp_si(res);
                return FLINT_MIN(acc, prec);
            }
            else if (arf_cmpabs_2exp_si(x, 1000) > 0)
            {
                if (fmpz_bits(ARF_EXPREF(x)) > 40)
                    return arb_lambertw_initial_asymp2(res, x, branch, prec);
                else
                    return arb_lambertw_initial_asymp1(res, x, branch, prec);
            }
            else
            {
                arf_set_d(res, d_lambertw(arf_get_d(x, ARF_RND_DOWN)));
                return 50;
            }
        }
        else
        {
            /* slightly smaller than the double exponent limit since
               d_lambertw_branch1 is unclever about underflowing */
            if (arf_cmpabs_2exp_si(x, -940) < 0)
            {
                if (fmpz_bits(ARF_EXPREF(x)) > 40)
                    return arb_lambertw_initial_asymp2(res, x, branch, prec);
                else
                    return arb_lambertw_initial_asymp1(res, x, branch, prec);
            }
            else
            {
                arf_set_d(res, d_lambertw_branch1(arf_get_d(x, ARF_RND_DOWN)));
                return 50;
            }
        }
    }
    else
    {
        /* Expand at -1/e */
        static const int coeffs[] = {-130636800,130636800,-43545600,19958400,
            -10402560,5813640,-3394560,2042589,-1256320};

        arb_t v;
        arf_t s;
        slong wp, k, acc;

        /* todo: could change precision dynamically depending on the
           closeness to -1/e */
        wp = 2 * prec + 20;

        arb_init(v);
        arf_init(s);

        arb_const_e(v, wp);
        arb_mul_arf(v, v, x, wp);
        arb_add_ui(v, v, 1, wp);
        arb_mul_2exp_si(v, v, 1);
        arb_sqrt(v, v, wp);
        if (branch)
            arb_neg(v, v);

        for (k = 8; k >= 0; k--)
        {
            arf_mul(s, s, arb_midref(v), wp, ARF_RND_DOWN);
            arf_add_si(s, s, coeffs[k], wp, ARF_RND_DOWN);
        }

        arf_div_si(s, s, -coeffs[0], wp, ARF_RND_DOWN);
        arf_set(res, s);

        arb_clear(v);
        arf_clear(s);

        /* Due to the arithmetic we should get no more than wp accurate bits */
        acc = wp;
        /* Truncation error is of order v^9 */
        if (!arf_is_special(arb_midref(v)))
            acc = FLINT_MIN(acc, 9 * (-ARF_EXP(arb_midref(v))));
        acc = FLINT_MAX(acc, 0);

        return acc;
    }
}

void
arb_lambertw(arb_t res, const arb_t x, int flags, slong prec)
{
    slong acc, wp, goal;
    slong ebits, ebits2;
    arb_t t, w;
    mag_t err;
    int branch;

    branch = flags & 1;

    if (!arb_is_finite(x))
    {
        arb_indeterminate(res);
        return;
    }

    if (branch == 1 && !arb_is_negative(x))
    {
        arb_indeterminate(res);
        return;
    }

    if (arb_is_zero(x))
    {
        arb_zero(res);
        return;
    }

    /* Quick estimate of log(x) and log(log(x)). */
    ebits = fmpz_bits(ARF_EXPREF(arb_midref(x)));
    ebits2 = FLINT_BIT_COUNT(ebits) + 2;

    /* Estimated accuracy goal. */
    goal = arb_rel_accuracy_bits(x);
    goal = FLINT_MAX(goal, 0);
    goal = FLINT_MIN(goal, prec);

    /* For huge x, we gain bits from the exponent. */
    if (branch == 0 && goal > 0 && arf_cmp_2exp_si(arb_midref(x), 10) > 0)
        goal = FLINT_MIN(goal + ebits - ebits2, prec);
    wp = goal + 4;

    /* Handle huge x directly. For x >= e,
        |W(x) - (log(x) - log(log(x)))| < 2/log(log(x))/log(x).
      (J. Inequal. Pure and Appl. Math., 9 (2) (2008), Art. 51).
      Note W(x) ~= log(x), so relative error is ~= log(log(x)) / log(x)^2. */
    if (branch == 0 && arf_cmp_2exp_si(arb_midref(x), 10) > 0
        && 2 * ebits - ebits2 > wp)
    {
        mag_t l, ll;
        arb_init(t);
        mag_init(l);
        mag_init(ll);

        arb_log(t, x, wp);
        arb_log(res, t, FLINT_MAX(wp - ebits + ebits2, 10));

        if (arb_is_positive(res))
        {
            arb_get_mag_lower(l, t);
            arb_get_mag(ll, res);
            arb_sub(res, t, res, prec);
            mag_div(l, ll, l);
            mag_mul_2exp_si(l, l, 1);
            arb_add_error_mag(res, l);
        }
        else
        {
            arb_indeterminate(res);
        }

        arb_clear(t);
        mag_clear(l);
        mag_clear(ll);
        return;
    }

    /* Handle tiny x directly. For k >= 2, |c_k| <= 4^k / 16. */
    if (branch == 0 && arf_cmpabs_2exp_si(arb_midref(x), -10) < 0
        && ebits > wp / 2)
    {
        mag_init(err);
        arb_get_mag(err, x);
        mag_mul_2exp_si(err, err, 2);

        if (ebits > wp)
        {
            arb_set_round(res, x, prec);
            mag_geom_series(err, err, 2);
        }
        else
        {
            arb_set(res, x);
            arb_submul(res, res, res, prec);
            mag_geom_series(err, err, 3);
        }

        mag_mul_2exp_si(err, err, -4);
        arb_add_error_mag(res, err);
        mag_clear(err);
        return;
    }

    arb_init(t);
    arb_init(w);
    mag_init(err);

    acc = arb_lambertw_initial(arb_midref(w), arb_midref(x), branch, wp);

    if (acc <= 2)
    {
        arb_indeterminate(w);
    }
    else if (acc >= wp)
    {
        arb_lambertw_bound_error(arb_radref(w), x, arb_midref(w),
            NULL, branch, FLINT_MAX(acc, 30));
    }
    else
    {
        slong k, padding, nextstep, maxstep, *steps;
        double rate, nearm1;

        steps = flint_malloc(sizeof(slong) * FLINT_BITS);

        /* Asymptotically, the Halley iteration triples the number of
           accurate bits. However, with a very large exponent we need
           some guard bits (due to evaluating the exponential?).
           This is heuristic. A better analysis should be possible. */
        rate = 2.0 + 1.0 / (1.0 + 0.01 * ebits);
        padding = 6 * ebits2;

        /* extra padding near -1/e */
        nearm1 = arf_get_d(arb_midref(w), ARF_RND_DOWN);
        if (fabs(nearm1 + 1.0) < 0.01)
        {
            arf_add_ui(arb_midref(t), arb_midref(w), 1, 30, ARF_RND_DOWN);

            if (arf_is_zero(arb_midref(t)))
                padding += prec;
            else
            {
                slong ee = -ARF_EXP(arb_midref(t));
                padding += FLINT_MIN(FLINT_MAX(2 * ee, 0), prec);
            }
        }

        maxstep = 0;
        steps[0] = wp;

        for (k = 1; k < FLINT_BITS; k++)
        {
            nextstep = steps[k - 1] / rate + padding;
            if (nextstep > acc)
            {
                steps[k] = nextstep;
                maxstep = k;
            }
            else
                break;
        }

        for (k = maxstep; k >= 0; k--)
        {
            arb_lambertw_halley_step(w, x, arb_midref(w),
                branch, (k == 0), steps[k] + 5);
        }

        flint_free(steps);
    }

    arb_set_round(res, w, prec);

    arb_clear(t);
    arb_clear(w);
    mag_clear(err);
}

