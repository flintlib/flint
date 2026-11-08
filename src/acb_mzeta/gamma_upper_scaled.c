/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"
#include "impl.h"

/* G(t) = t^-a Gamma(a, t) = int_1^oo u^(a-1) e^(-t u) du for real t > 0.

   The lattice sums need G at many points with one a, to a relative
   accuracy chosen by the caller, so the choice of method and of the
   working precision is made here rather than delegated to
   acb_hypgeom_gamma_upper (whose tuning targets other uses):

   - asymptotic series: F(t) = e^t G(t) = int_0^oo (1+v)^(a-1) e^(-tv) dv
     = sum_{k<N} (a-1)...(a-k) / t^(k+1) + R_N with
     |R_N| <= |(a-1)...(a-N)| / t^(N+1) for N >= Re(a) - 1 (Taylor remainder
     of (1+v)^(a-1) in integral form, |(1 + theta v)^(a-1-N)| <= 1);
     usable when the smallest term is below the tolerance (t of order
     prec log 2 or larger);
   - convergent series: G(t) = Gamma(a) t^-a - e^-t sum_k t^k / (a)_{k+1},
     whose cancellation (about t / log 2 bits, more near the poles of
     Gamma) is compensated by raising the working precision by an
     estimate of the loss computed in double arithmetic;
   - a = -n a nonpositive integer: G(t) = t^n Gamma(-n, t) =
     t^n (-1)^n / n! (E1(t) - e^-t sum_{k<n} (-1)^k k! / t^(k+1)) with
     E1(t) = -gamma - log t + e^-t sum_{k>=1} H_k t^k / k! (positive terms).

   A ball t is handled by evaluating at its midpoint and adding
   rad(t) sup |G'| with |G'(s)| = |G_{a+1}(s)| <= G_{Re(a)+1}(lo t), which
   avoids the amplification of the input radius (by about e^t) that the
   series would cause. */

#define LOG2E 1.4426950408889634

int
_acb_mzeta_gamma_upper_scaled_asymp(acb_t res, const acb_t a, const arb_t t, slong prec)
{
    double td, ar, ai, lrel, ltd;
    arf_t aru;
    slong k, N = -1, wp;
    acb_t s, term, u;
    arb_t tinv;
    mag_t err;

    td = arf_get_d(arb_midref(t), ARF_RND_DOWN);
    if (!(td > 0.0) || !arb_is_exact(t))
        return 0;

    ar = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    ai = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    arf_init(aru);
    arb_get_ubound_arf(aru, acb_realref(a), 30);
    ltd = log2(td);

    /* log2 of |term_k| t = |(a-1)...(a-k)| / t^k */
    lrel = 0.0;
    for (k = 1; k < 1000000; k++)
    {
        double ak = hypot(ar - k, ai);
        lrel += log2(FLINT_MAX(ak, 1e-300)) - ltd;
        if (lrel < -(double) prec - 4 && arf_cmp_si(aru, k + 1) <= 0)
        {
            N = k;
            break;
        }
        /* past the smallest term */
        if (ak > td && k > ar + 1)
            break;
    }
    arf_clear(aru);

    if (N == -1)
        return 0;

    wp = prec + 10 + FLINT_BIT_COUNT(N);
    acb_init(s);
    acb_init(term);
    acb_init(u);
    arb_init(tinv);
    mag_init(err);

    arb_inv(tinv, t, wp);
    acb_set_arb(term, tinv);
    for (k = 0; k < N; k++)
    {
        acb_add(s, s, term, wp);
        /* term_{k+1} = term_k (a - k - 1) / t */
        acb_sub_ui(u, a, k + 1, wp);
        acb_mul(term, term, u, wp);
        acb_mul_arb(term, term, tinv, wp);
    }
    acb_get_mag(err, term);
    acb_add_error_mag(s, err);
    acb_swap(res, s);

    acb_clear(s);
    acb_clear(term);
    acb_clear(u);
    arb_clear(tinv);
    mag_clear(err);
    return 1;
}

/* G for a = -n, n >= 0 an integer */
static void
gamma_upper_scaled_negint(acb_t res, ulong n, const arb_t t, slong prec)
{
    double td = arf_get_d(arb_midref(t), ARF_RND_DOWN);
    slong extra, wp, k;
    arb_t s, p, H, v, lg, e;
    mag_t err;

    /* cancellation in E1 (about t log2 e bits) and in the bracket */
    if (!(td > 0.0 && td < 1e12))
    {
        acb_indeterminate(res);
        return;
    }
    extra = (slong) (td * LOG2E + log2(fabs(log(td)) + 2.0)
                     + n * (fabs(log2(td)) + log2(n + 1.0))) + 20;
    wp = prec + FLINT_MAX(extra, 0);

    arb_init(s); arb_init(p); arb_init(H); arb_init(v); arb_init(lg); arb_init(e);
    mag_init(err);

    /* sum_{k>=1} H_k t^k / k! */
    arb_one(p);
    arb_zero(H);
    for (k = 1; ; k++)
    {
        arb_mul(p, p, t, wp);
        arb_div_ui(p, p, k, wp);
        arb_set_ui(v, 1);
        arb_div_ui(v, v, k, wp);
        arb_add(H, H, v, wp);
        arb_mul(v, p, H, wp);
        arb_add(s, s, v, wp);
        /* for k >= 2t the ratios of p are <= 1/2 and H grows by at most
           1/(k+1) per step: the tail is at most 2 p_k (H_k + 1) */
        /* for j > k >= 2t - 1, p_j / p_{j-1} = t / j <= 1/2 and
           H_j <= H_k + (j - k), so the tail is at most p_k (H_k + 2) */
        if (k + 1 >= 2 * td + 1 && arf_cmpabs_2exp_si(arb_midref(v), arf_abs_bound_lt_2exp_si(arb_midref(s)) - wp - 4) < 0)
        {
            arb_mul_2exp_si(lg, t, 1);
            arb_sub_ui(lg, lg, k + 1, wp);
            if (arb_is_nonpositive(lg))
            {
                arb_add_ui(v, H, 2, wp);
                arb_mul(v, v, p, wp);
                arb_get_mag(err, v);
                arb_add_error_mag(s, err);
                break;
            }
        }
        if (k > 100 * (td + 10) + 10 * wp)
        {
            arb_indeterminate(s);
            break;
        }
    }

    /* E1 = -gamma - log t + e^-t s */
    arb_neg(e, t);
    arb_exp(e, e, wp);
    arb_mul(s, s, e, wp);
    arb_const_euler(v, wp);
    arb_sub(s, s, v, wp);
    arb_log(lg, t, wp);
    arb_sub(s, s, lg, wp);

    /* - e^-t sum_{k<n} (-1)^k k! / t^(k+1) */
    if (n > 0)
    {
        arb_inv(p, t, wp);
        arb_set(v, p);          /* 0! / t */
        arb_zero(H);
        for (k = 0; k < (slong) n; k++)
        {
            if (k % 2 == 0)
                arb_add(H, H, v, wp);
            else
                arb_sub(H, H, v, wp);
            arb_mul_ui(v, v, k + 1, wp);
            arb_mul(v, v, p, wp);
        }
        arb_mul(H, H, e, wp);
        arb_sub(s, s, H, wp);
    }

    /* t^n (-1)^n / n! */
    arb_pow_ui(v, t, n, wp);
    arb_mul(s, s, v, wp);
    arb_fac_ui(v, n, wp);
    arb_div(s, s, v, wp);
    if (n % 2)
        arb_neg(s, s);

    acb_set_arb(res, s);
    acb_set_round(res, res, prec);

    arb_clear(s); arb_clear(p); arb_clear(H); arb_clear(v); arb_clear(lg); arb_clear(e);
    mag_clear(err);
}

/* S = sum_k t^k / (a)_{k+1} = (1/a) sum_k t^k / (a+1)_k for exact dyadic
   a > -1, a != 0, by _arb_hypgeom_gamma_lower_sum_rs_1 with a + 1 = p/q;
   returns 0 if a does not have this form.  The sum stops at N with
   a + 1 + N >= 2t, where the ratios t / (a + 1 + k) are at most 1/2, so the
   tail is at most twice the term of index N. */
static int
series_real_rs(arb_t S, const arb_t a, double ar, const arb_t t, double td, slong wp)
{
    fmpz_t m, e;
    ulong p, q;
    slong N, k;
    double lterm, lsum, l2t;
    arb_t a1, u, v;
    mag_t err;
    int ok;

    if (!arb_is_exact(a) || arf_is_zero(arb_midref(a)) || !(ar > -1.0) || !(td < 1e6))
        return 0;

    fmpz_init(m);
    fmpz_init(e);
    arf_get_fmpz_2exp(m, e, arb_midref(a));
    ok = fmpz_cmp_si(e, 0) <= 0 && fmpz_cmp_si(e, -(FLINT_BITS / 2)) >= 0;
    if (ok)
    {
        q = UWORD(1) << (-fmpz_get_si(e));
        fmpz_add_ui(m, m, q);           /* (a + 1) q */
        ok = fmpz_sgn(m) > 0 && fmpz_bits(m) < FLINT_BITS / 2;
        p = ok ? fmpz_get_ui(m) : 0;
    }
    fmpz_clear(m);
    fmpz_clear(e);
    if (!ok)
        return 0;

    /* N: past 2t - a - 1, and the term below 2^-wp times the sum (whose
       largest term is a lower bound) */
    l2t = log2(td);
    lterm = 0.0;
    lsum = 0.0;
    for (N = 1; ; N++)
    {
        lterm += l2t - log2(ar + N);
        lsum = FLINT_MAX(lsum, lterm);
        if (ar + 1 + N >= 2 * td + 1 && lterm < lsum - wp - 2)
            break;
        if (N > 100 * (td + 10) + 10 * wp)
            return 0;
    }
    if (!((double) q * (N + 1) + p < ldexp(1.0, FLINT_BITS - 2)))
        return 0;

    arb_init(a1);
    arb_init(u);
    arb_init(v);
    mag_init(err);

    _arb_hypgeom_gamma_lower_sum_rs_1(S, p, q, t, N, wp);

    /* the tail: 2 t^N / (a+1)_N (with a + 1 + N >= 2t, checked) */
    arb_add_ui(a1, a, 1, 30);
    arb_add_ui(v, a1, N, 30);
    arb_submul_ui(v, t, 2, 30);
    ok = arb_is_nonnegative(v);
    if (ok)
    {
        arb_pow_ui(u, t, N, 30);
        arb_hypgeom_rising_ui(v, a1, N, 30);
        arb_div(u, u, v, 30);
        arb_get_mag(err, u);
        mag_mul_2exp_si(err, err, 1);
        arb_add_error_mag(S, err);
        arb_div(S, S, a, wp);
    }

    arb_clear(a1);
    arb_clear(u);
    arb_clear(v);
    mag_clear(err);
    (void) k;
    return ok;
}

/* the convergent series for real a, in real arithmetic (wp from the
   estimates of gamma_upper_scaled_series) */
static void
series_real(acb_t res, const arb_t a, double ar, const arb_t t, double td, const arb_t ga,
    slong prec, slong wp)
{
    arb_t s, u, v, g, lt;
    slong k;

    arb_init(s); arb_init(u); arb_init(v); arb_init(g); arb_init(lt);

    /* a + 1 = p/q > 0 with q a power of two: S = (1/a) sum_k t^k / (a+1)_k
       by rectangular splitting (no division per term) */
    if (series_real_rs(s, a, ar, t, td, wp))
        goto assemble;

    /* S = sum_k t^k / (a)_{k+1}: u_0 = 1/a, u_{k+1} = u_k t / (a + k + 1) */
    arb_inv(u, a, wp);
    arb_zero(s);
    for (k = 0; ; k++)
    {
        arb_add(s, s, u, wp);
        arb_add_ui(v, a, k + 1, wp);
        arb_mul(u, u, t, wp);
        arb_div(u, u, v, wp);
        /* once a + k + 2 >= 2t, the remaining ratios are <= 1/2 */
        if (ar + k + 2 >= 2 * td + 1 && arb_is_finite(s) && arb_is_finite(u))
        {
            mag_t mu, ms;
            int small;
            mag_init(mu);
            mag_init(ms);
            arb_get_mag(mu, u);
            arb_get_mag_lower(ms, s);
            mag_mul_2exp_si(ms, ms, -wp);
            small = mag_cmp(mu, ms) <= 0;
            mag_clear(ms);
            if (small)
            {
                /* rigorous: a + k + 2 >= 2t, so |t / (a + j + 1)| <= 1/2
                   for j >= k + 1, and the tail is at most 2 |u_{k+1}| */
                arb_add_ui(lt, a, k + 2, wp);
                arb_submul_ui(lt, t, 2, wp);
                if (arb_is_nonnegative(lt))
                {
                    mag_mul_2exp_si(mu, mu, 1);
                    arb_add_error_mag(s, mu);
                    mag_clear(mu);
                    break;
                }
            }
            mag_clear(mu);
        }
        if (k > 100 * (td + 10) + 10 * wp)
        {
            arb_indeterminate(s);
            break;
        }
    }

assemble:
    /* e^-t S */
    arb_neg(lt, t);
    arb_exp(lt, lt, wp);
    arb_mul(s, s, lt, wp);

    /* Gamma(a) t^-a */
    if (ga != NULL && arb_rel_accuracy_bits(ga) >= wp)
        arb_set(g, ga);
    else
        arb_gamma(g, a, wp);
    arb_log(lt, t, wp);
    arb_mul(v, a, lt, wp);
    arb_neg(v, v);
    arb_exp(v, v, wp);
    arb_mul(g, g, v, wp);

    arb_sub(acb_realref(res), g, s, wp);
    arb_set_round(acb_realref(res), acb_realref(res), prec);
    arb_zero(acb_imagref(res));

    arb_clear(s); arb_clear(u); arb_clear(v); arb_clear(g); arb_clear(lt);
}

/* G at an exact t > 0 by the convergent series */
static void
gamma_upper_scaled_series(acb_t res, const acb_t a, const arb_t t, const acb_t ga, slong prec)
{
    double td, ar, ai, lmax, lterm, lg, lgam, lres;
    slong k, extra, wp;
    acb_t s, u, v, g;
    arb_t lt;
    mag_t err;

    td = arf_get_d(arb_midref(t), ARF_RND_DOWN);
    ar = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    ai = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);

    acb_init(s); acb_init(u); acb_init(v); acb_init(g);
    arb_init(lt);
    mag_init(err);

    /* estimates (log2): the largest term of e^-t S, |Gamma(a) t^-a|, and a
       lower estimate of |G| ~ e^-t / (t + |1 - a| + 1); the terms
       t^k / |(a)_{k+1}| decrease once |a + k + 1| > t, so the loop stops
       there */
    {
        double l2t = log2(td);
        lterm = -log2(FLINT_MAX(hypot(ar, ai), 1e-300));
        lmax = lterm;
        for (k = 1; k < 4 * td + 100; k++)
        {
            lterm += l2t - log2(FLINT_MAX(hypot(ar + k, ai), 1e-300));
            lmax = FLINT_MAX(lmax, lterm);
            if (ar + k + 1 > td + 1)
                break;
        }
    }
    lmax -= td * LOG2E;
    if (ga != NULL && acb_is_finite(ga) && !acb_contains_zero(ga))
    {
        /* log2 |Gamma(a)| from the given value */
        mag_t m;
        mag_init(m);
        acb_get_mag(m, ga);
        lg = mag_get_d_log2_approx(m);
        mag_clear(m);
    }
    else
    {
        acb_lgamma(g, a, 30);
        lg = arf_get_d(arb_midref(acb_realref(g)), ARF_RND_NEAR) * LOG2E;
        if (!arf_is_finite(arb_midref(acb_realref(g))))
            lg = 0.0;
    }
    lgam = lg - ar * log2(td);
    lres = -td * LOG2E - log2(td + hypot(1.0 - ar, ai) + 1.0);
    if (!(td > 0.0 && td < 1e12) || !(fabs(lmax) < 1e15 && fabs(lgam) < 1e15))
    {
        acb_indeterminate(res);
        goto cleanup;
    }
    extra = (slong) (FLINT_MAX(lmax, lgam) - lres) + 12;
    wp = prec + FLINT_MAX(extra, 0) + 8;

    if (arb_is_zero(acb_imagref(a)) && (ga == NULL || arb_is_zero(acb_imagref(ga))))
    {
        series_real(res, acb_realref(a), ar, t, td, ga == NULL ? NULL : acb_realref(ga), prec, wp);
        goto cleanup;
    }

    /* S = sum_k t^k / (a)_{k+1}: u_0 = 1/a, u_{k+1} = u_k t / (a + k + 1) */
    acb_inv(u, a, wp);
    acb_zero(s);
    for (k = 0; ; k++)
    {
        acb_add(s, s, u, wp);
        acb_add_ui(v, a, k + 1, wp);
        acb_div(u, u, v, wp);
        acb_mul_arb(u, u, t, wp);
        /* once Re(a) + k + 2 >= 2t, the remaining ratios are <= 1/2 */
        if (ar + k + 2 >= 2 * td + 1 && acb_is_finite(s) && acb_is_finite(u))
        {
            mag_t mu, ms;
            int small;
            mag_init(mu);
            mag_init(ms);
            acb_get_mag(mu, u);
            acb_get_mag_lower(ms, s);
            mag_mul_2exp_si(ms, ms, -wp);
            small = mag_cmp(mu, ms) <= 0;
            mag_clear(ms);
            if (small)
            {
                /* rigorous: Re(a) + k + 2 >= 2t, so |t / (a + j + 1)| <= 1/2
                   for j >= k + 1, and the tail is at most 2 |u_{k+1}| */
                arb_add_ui(lt, acb_realref(a), k + 2, wp);
                arb_submul_ui(lt, t, 2, wp);
                if (arb_is_nonnegative(lt))
                {
                    mag_mul_2exp_si(mu, mu, 1);
                    acb_add_error_mag(s, mu);
                    mag_clear(mu);
                    break;
                }
            }
            mag_clear(mu);
        }
        if (k > 100 * (td + 10) + 10 * wp)
        {
            acb_indeterminate(s);
            break;
        }
    }

    /* e^-t S */
    arb_neg(lt, t);
    arb_exp(lt, lt, wp);
    acb_mul_arb(s, s, lt, wp);

    /* Gamma(a) t^-a */
    if (ga != NULL && acb_rel_accuracy_bits(ga) >= wp)
        acb_set(g, ga);
    else
        acb_gamma(g, a, wp);
    arb_log(lt, t, wp);
    acb_mul_arb(v, a, lt, wp);
    acb_neg(v, v);
    acb_exp(v, v, wp);
    acb_mul(g, g, v, wp);

    acb_sub(res, g, s, wp);
    acb_set_round(res, res, prec);

cleanup:
    acb_clear(s); acb_clear(u); acb_clear(v); acb_clear(g);
    arb_clear(lt);
    mag_clear(err);
}

static void
gamma_upper_scaled_exact(acb_t res, const acb_t a, const arb_t t, const acb_t ga, slong prec)
{
    if (_acb_mzeta_gamma_upper_scaled_asymp(res, a, t, prec))
    {
        arb_t e;
        arb_init(e);
        arb_neg(e, t);
        arb_exp(e, e, prec + 10);
        acb_mul_arb(res, res, e, prec + 10);
        acb_set_round(res, res, prec);
        arb_clear(e);
    }
    else if (acb_is_int(a) && arb_is_nonpositive(acb_realref(a)))
    {
        fmpz_t n;
        fmpz_init(n);
        arf_get_fmpz(n, arb_midref(acb_realref(a)), ARF_RND_NEAR);
        fmpz_neg(n, n);
        if (fmpz_cmp_ui(n, 1000000) < 0)
            gamma_upper_scaled_negint(res, fmpz_get_ui(n), t, prec);
        else
            acb_indeterminate(res);
        fmpz_clear(n);
    }
    else
    {
        gamma_upper_scaled_series(res, a, t, ga, prec);
    }
}

void
_acb_mzeta_gamma_upper_scaled(acb_t res, const acb_t a, const arb_t t, const acb_t ga, slong prec)
{
    if (!arb_is_positive(t) || !acb_is_finite(a))
    {
        acb_indeterminate(res);
        return;
    }

    if (arb_is_exact(t))
    {
        gamma_upper_scaled_exact(res, a, t, ga, prec);
    }
    else
    {
        arb_t tm, lo;
        acb_t b;
        mag_t err;

        arb_init(tm);
        arb_init(lo);
        acb_init(b);
        mag_init(err);

        arb_set_arf(tm, arb_midref(t));
        gamma_upper_scaled_exact(res, a, tm, ga, prec);

        /* rad(t) G_{Re(a) + 1}(lo), G_{a'} decreasing on (0, oo) for real a' */
        arb_get_lbound_arf(arb_midref(lo), t, 30);
        arb_set_arf(lo, arb_midref(lo));
        arb_add_ui(acb_realref(b), acb_realref(a), 1, 30);
        arb_get_ubound_arf(arb_midref(acb_realref(b)), acb_realref(b), 30);
        mag_zero(arb_radref(acb_realref(b)));
        if (arf_sgn(arb_midref(lo)) > 0)
        {
            acb_t g;
            acb_init(g);
            gamma_upper_scaled_exact(g, b, lo, NULL, 30);
            acb_get_mag(err, g);
            mag_mul(err, err, arb_radref(t));
            acb_add_error_mag(res, err);
            acb_clear(g);
        }
        else
            acb_indeterminate(res);

        arb_clear(tm);
        arb_clear(lo);
        acb_clear(b);
        mag_clear(err);
    }
}
