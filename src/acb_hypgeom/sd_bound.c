/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Rigorous bounds for the integrand exp(phi(w) - c),

        phi(w) = alpha e^w + beta e^(-w) - nu w,

    with beta = -alpha (sym) or beta = 0, for all alpha and nu in the
    given balls.  Everything that enters an error bound or a bound
    returned to the user is computed here, in arb and mag arithmetic.
    The contours themselves, the choice of the quadrature rules and of
    the parameters of the expansions, and all tolerances are chosen
    heuristically in double precision elsewhere (sd_plan.c,
    sd_integrate.c, sd_expand.c); since the integrand is entire, they
    affect only the size of the bounds computed here, never their
    validity.  The only data taken over from the heuristic code are
    exact: the vertices of the polygons and the reference values c
    (doubles), and the kinds of the ends of the polygons (the valley
    which a tail is to end in is verified here).

    The three weights are g_0 = 1, g_1 = e^x + [beta != 0] e^(-x) and
    g_2 = |w| (x = Re w); the weighted integrals bound the derivatives
    of the integral with respect to alpha and nu.
*/

#include "acb.h"
#include "acb_hypgeom/sd.h"

/* The working precision.  The values of phi are of magnitude up to
   2^44 e^|x| or so and are needed with absolute errors well below 1;
   the cost in arb arithmetic hardly depends on the precision here. */
#define SDB_PREC 128

/* phi, phi', phi'' and phi''' = phi' + nu at w (the last computed
   without the cancellation; d3 can be NULL) */
void
_sd_bound_phi(acb_t phi, acb_t d1, acb_t d2, acb_t d3, const acb_t alpha, int sym,
    const acb_t nu, const acb_t w)
{
    acb_t a, b;

    acb_init(a);
    acb_init(b);

    acb_exp(a, w, SDB_PREC);
    if (sym)
    {
        /* b = beta e^-w = -alpha e^-w */
        acb_div(b, alpha, a, SDB_PREC);
        acb_neg(b, b);
    }
    acb_mul(a, a, alpha, SDB_PREC);

    acb_mul(phi, nu, w, SDB_PREC);
    acb_sub(d1, a, b, SDB_PREC);
    if (d3 != NULL)
        acb_set(d3, d1);
    acb_sub(d1, d1, nu, SDB_PREC);
    acb_add(d2, a, b, SDB_PREC);
    acb_sub(phi, d2, phi, SDB_PREC);

    acb_clear(a);
    acb_clear(b);
}

/* res = an upper bound for exp(x) */
static void
sdb_exp_upper(mag_t res, const arb_t x)
{
    arb_t t;
    arb_init(t);
    arb_exp(t, x, SDB_PREC);
    arb_get_mag(res, t);
    arb_clear(t);
}

/* upper bounds for cosh R - 1 and sinh R - R */
static void
sdb_cosh_sinh_rem(mag_t chm1, mag_t shmr, const mag_t R)
{
    if (mag_cmp_2exp_si(R, -4) < 0)
    {
        /* cosh R - 1 <= (R^2 / 2) cosh R, sinh R - R <= (R^3 / 6) cosh R
           (termwise) */
        mag_t c;
        mag_init(c);
        mag_cosh(c, R);
        mag_mul(chm1, R, R);
        mag_mul(shmr, chm1, R);
        mag_mul(chm1, chm1, c);
        mag_mul_2exp_si(chm1, chm1, -1);
        mag_mul(shmr, shmr, c);
        mag_div_ui(shmr, shmr, 6);
        mag_clear(c);
    }
    else
    {
        arb_t r, s, c;
        arb_init(r); arb_init(s); arb_init(c);
        arf_set_mag(arb_midref(r), R);
        arb_sinh_cosh(s, c, r, SDB_PREC);
        arb_sub_ui(c, c, 1, SDB_PREC);
        arb_sub(s, s, r, SDB_PREC);
        arb_get_mag(chm1, c);
        arb_get_mag(shmr, s);
        arb_clear(r); arb_clear(s); arb_clear(c);
    }
}

/* Since alpha e^(m+u) + beta e^(-m-u) = phi''(m) cosh u + phi'''(m) sinh u,

       phi(m + u) - phi(m) = phi'(m) u + phi''(m) (cosh u - 1)
                             + phi'''(m) (sinh u - u),

   so that for |u| <= R

       Re(phi(m + u) - c) <= Re(phi(m) - c) + |phi'(m)| R
           + |phi''(m)| (cosh R - 1) + |phi'''(m)| (sinh R - R). */
void
_sd_bound_disk(mag_t M, const acb_t alpha, int sym, const acb_t nu,
    const acb_t m, const mag_t R, const acb_t c)
{
    acb_t phi, d1, d2, d3;
    mag_t t, u, chm1, shmr;
    arb_t e;

    acb_init(phi); acb_init(d1); acb_init(d2); acb_init(d3);
    mag_init(t); mag_init(u); mag_init(chm1); mag_init(shmr);
    arb_init(e);

    _sd_bound_phi(phi, d1, d2, d3, alpha, sym, nu, m);
    sdb_cosh_sinh_rem(chm1, shmr, R);

    acb_get_mag(t, d1);
    mag_mul(t, t, R);
    acb_get_mag(u, d2);
    mag_addmul(t, u, chm1);
    acb_get_mag(u, d3);
    mag_addmul(t, u, shmr);

    arb_sub(e, acb_realref(phi), acb_realref(c), SDB_PREC);
    arb_add_error_mag(e, t);
    /* (the upper bound of the ball is what matters) */
    sdb_exp_upper(M, e);

    acb_clear(phi); acb_clear(d1); acb_clear(d2); acb_clear(d3);
    mag_clear(t); mag_clear(u); mag_clear(chm1); mag_clear(shmr);
    arb_clear(e);
}

/* upper bounds for the suprema of the weights g_1 and g_2 on the
   segment [a, b] */
static void
sdb_weights(mag_t g1, mag_t g2, int sym, const acb_t a, const acb_t b)
{
    arb_t t;
    mag_t m;

    arb_init(t);
    mag_init(m);

    /* e^x <= e^max(Re a, Re b), e^-x <= e^max(-Re a, -Re b) */
    arb_max(t, acb_realref(a), acb_realref(b), SDB_PREC);
    sdb_exp_upper(g1, t);
    if (sym)
    {
        arb_min(t, acb_realref(a), acb_realref(b), SDB_PREC);
        arb_neg(t, t);
        sdb_exp_upper(m, t);
        mag_add(g1, g1, m);
    }

    /* |w| is convex */
    acb_get_mag(g2, a);
    acb_get_mag(m, b);
    mag_max(g2, g2, m);

    arb_clear(t);
    mag_clear(m);
}

/* With w = wm + h t, -1 <= t <= 1 (wm = (a + b) / 2, h = (b - a) / 2),

       Re(phi(w) - c) <= A + B t + C t^2,
       A = Re(phi(wm) - c) + |h|^3 S / 6,  B = Re(phi'(wm) h),
       C = Re(phi''(wm) h^2) / 2,

   where S >= |alpha| e^x + |beta| e^-x >= |phi'''| on the segment.
   With b >= |B| and c >= C, the integral of exp(A + B t + C t^2) over
   [-1, 1] is at most 2 exp(A + max_{0 <= s <= 1} (b s + c s^2)) and,
   if c < 0, at most sqrt(pi / -c) exp(A + b^2 / (-4 c)). */
void
_sd_bound_segment(mag_ptr res, int nw, const acb_t alpha, int sym, const acb_t nu,
    const acb_t a, const acb_t b, const acb_t c)
{
    acb_t wm, h, phi, d1, d2, t;
    arb_t A, bb, cc, E, u, v;
    mag_t S, ah, g1, g2, val, m;

    acb_init(wm); acb_init(h); acb_init(phi); acb_init(d1); acb_init(d2); acb_init(t);
    arb_init(A); arb_init(bb); arb_init(cc); arb_init(E); arb_init(u); arb_init(v);
    mag_init(S); mag_init(ah); mag_init(g1); mag_init(g2); mag_init(val); mag_init(m);

    acb_add(wm, a, b, SDB_PREC);
    acb_mul_2exp_si(wm, wm, -1);
    acb_sub(h, b, a, SDB_PREC);
    acb_mul_2exp_si(h, h, -1);
    acb_get_mag(ah, h);

    _sd_bound_phi(phi, d1, d2, NULL, alpha, sym, nu, wm);
    sdb_weights(g1, g2, sym, a, b);

    /* S = |alpha| (e^xhi + [sym] e^-xlo), with g1 for the second factor */
    acb_get_mag(S, alpha);
    mag_mul(S, S, g1);

    /* A */
    mag_mul(m, ah, ah);
    mag_mul(m, m, ah);
    mag_mul(m, m, S);
    mag_div_ui(m, m, 6);
    arb_sub(A, acb_realref(phi), acb_realref(c), SDB_PREC);
    arb_add_error_mag(A, m);
    arb_get_ubound_arf(arb_midref(A), A, SDB_PREC);
    mag_zero(arb_radref(A));

    /* b >= |B| */
    acb_mul(t, d1, h, SDB_PREC);
    arb_get_mag(m, acb_realref(t));
    arf_set_mag(arb_midref(bb), m);

    /* c >= C */
    acb_mul(t, h, h, SDB_PREC);
    acb_mul(t, t, d2, SDB_PREC);
    arb_mul_2exp_si(cc, acb_realref(t), -1);
    arb_get_ubound_arf(arb_midref(cc), cc, SDB_PREC);
    mag_zero(arb_radref(cc));

    if (!arb_is_finite(A) || !arb_is_finite(bb) || !arb_is_finite(cc))
    {
        mag_inf(val);
    }
    else
    {
        int neg = arf_sgn(arb_midref(cc)) < 0;

        /* u = b^2 / (-4 c) if c < 0 */
        if (neg)
        {
            arb_mul(u, bb, bb, SDB_PREC);
            arb_div(u, u, cc, SDB_PREC);
            arb_mul_2exp_si(u, u, -2);
            arb_neg(u, u);
        }

        /* The maximum of b s + c s^2 on [0, 1] is b + c if c >= 0 or
           if c < 0 and b + 2 c >= 0; otherwise it is the maximum over
           all s, u = b^2 / (-4 c) (which is a bound in any case, and
           is used if the sign of b + 2 c is not decided). */
        if (neg)
        {
            arb_mul_2exp_si(v, cc, 1);
            arb_add(v, v, bb, SDB_PREC);
            if (arb_is_nonnegative(v))
                arb_add(E, bb, cc, SDB_PREC);
            else
                arb_set(E, u);
        }
        else
            arb_add(E, bb, cc, SDB_PREC);
        arb_add(E, E, A, SDB_PREC);
        sdb_exp_upper(val, E);
        mag_mul_2exp_si(val, val, 1);

        if (neg)
        {
            arb_const_pi(v, SDB_PREC);
            arb_div(v, v, cc, SDB_PREC);
            arb_neg(v, v);
            arb_sqrt(v, v, SDB_PREC);
            arb_add(E, A, u, SDB_PREC);
            arb_exp(E, E, SDB_PREC);
            arb_mul(E, E, v, SDB_PREC);
            arb_get_mag(m, E);
            mag_min(val, val, m);
        }
    }

    /* |dw| = |h| dt */
    mag_mul(val, val, ah);

    if (nw > 0)
        mag_set(res, val);
    if (nw > 1)
        mag_mul(res + 1, val, g1);
    if (nw > 2)
        mag_mul(res + 2, val, g2);

    acb_clear(wm); acb_clear(h); acb_clear(phi); acb_clear(d1); acb_clear(d2); acb_clear(t);
    arb_clear(A); arb_clear(bb); arb_clear(cc); arb_clear(E); arb_clear(u); arb_clear(v);
    mag_clear(S); mag_clear(ah); mag_clear(g1); mag_clear(g2); mag_clear(val); mag_clear(m);
}

/* res >= int_1^inf e^(-K (u - 1)) u^b du for K > 0, given a lower
   bound for K and an upper bound for b (exact numbers as balls):
   1 / K if b <= 0; 1 / (K - b) if 0 < b < K, since u^b <= e^(b (u - 1));
   and otherwise, since u^b e^(-K (u - 1) / 2) is largest at
   u = 2 b / K >= 1,

       (2 / K) (2 b / K)^b e^(K / 2 - b). */
static void
sdb_gamma_tail(mag_t res, const arb_t K, const arb_t b)
{
    arb_t t, u;

    arb_init(t);
    arb_init(u);

    arb_sub(u, K, b, SDB_PREC);

    if (!arb_is_positive(b))
    {
        /* (b is an exact number: b <= 0) */
        arb_inv(t, K, SDB_PREC);
    }
    else if (arb_is_positive(u))
    {
        arb_inv(t, u, SDB_PREC);
    }
    else
    {
        /* K <= b, or undecided (then K - b is tiny, and the bound
           holds as 2 b / K >= 1 up to the precision: make sure) */
        arb_mul_2exp_si(t, b, 1);
        arb_sub(u, t, K, SDB_PREC);
        if (!arb_is_nonnegative(u))
        {
            arb_pos_inf(t);
        }
        else
        {
            arb_div(t, t, K, SDB_PREC);         /* 2 b / K */
            arb_log(t, t, SDB_PREC);
            arb_mul(t, t, b, SDB_PREC);
            arb_mul_2exp_si(u, K, -1);
            arb_add(t, t, u, SDB_PREC);
            arb_sub(t, t, b, SDB_PREC);
            arb_exp(t, t, SDB_PREC);
            arb_div(t, t, K, SDB_PREC);
            arb_mul_2exp_si(t, t, 1);
        }
    }

    arb_get_mag(res, t);

    arb_clear(t);
    arb_clear(u);
}

/* The ray w = we + s x', x' >= 0 (s = 1: to the valley SD_PLUS, s = -1:
   to SD_MINUS).  With u = e^x' >= 1,

       Re phi(w) = P u + Q / u - Re(nu we) - s Re(nu) log u,

   where P = Re(alpha e^we), Q = Re(beta e^-we) for s = 1 and the other
   way round for s = -1.  If P <= -K < 0 for all alpha in the ball,

       Re(phi(w) - c) <= G - K (u - 1) + a log u,
       G >= Re(phi(we) - c),   a = max(-Q, 0) + max(-s Re nu, 0)

   (using Q / u - Q <= max(-Q, 0) log u), and |dw| = du / u.  The
   weights are at most g_1(we) u and (|we| + 1) u on the ray.  The ray
   stays in the valley with index k: the integer nearest to
   (Im we + arg(alpha) - pi) / (2 pi) for s = 1, and to
   (Im we - arg(-beta)) / (2 pi) for s = -1, which is unique if P < 0. */
static int
sdb_tail_exp(mag_ptr res, int nw, int * kout, const acb_t alpha, int sym, const acb_t nu,
    cd we, const acb_t c, int s)
{
    acb_t w, e, t;
    arb_t P, Q, K, a, G, u, v;
    mag_t eg, m, g;
    int ok = 0;

    if (s == -1 && !sym)
        return 0;

    acb_init(w); acb_init(e); acb_init(t);
    arb_init(P); arb_init(Q); arb_init(K); arb_init(a); arb_init(G); arb_init(u); arb_init(v);
    mag_init(eg); mag_init(m); mag_init(g);

    arb_set_d(acb_realref(w), we.re);
    arb_set_d(acb_imagref(w), we.im);

    /* P and Q */
    acb_exp(e, w, SDB_PREC);
    acb_mul(t, alpha, e, SDB_PREC);
    arb_set(P, acb_realref(t));
    arb_zero(Q);
    if (sym)
    {
        acb_div(t, alpha, e, SDB_PREC);
        arb_neg(Q, acb_realref(t));
    }
    if (s == -1)
        arb_swap(P, Q);

    /* K = a lower bound for -P */
    arb_neg(K, P);
    arb_get_lbound_arf(arb_midref(K), K, SDB_PREC);
    mag_zero(arb_radref(K));
    if (!arb_is_finite(K) || !arb_is_positive(K))
        goto cleanup;

    /* the valley */
    acb_arg(u, alpha, SDB_PREC);
    arb_const_pi(v, SDB_PREC);
    if (s == 1)
    {
        arb_add(u, u, acb_imagref(w), SDB_PREC);
        arb_sub(u, u, v, SDB_PREC);
    }
    else
        arb_sub(u, acb_imagref(w), u, SDB_PREC);
    arb_div(u, u, v, SDB_PREC);
    arb_mul_2exp_si(u, u, -1);
    if (!arb_is_finite(u) || arf_cmpabs_2exp_si(arb_midref(u), 30) >= 0)
        goto cleanup;
    *kout = (int) floor(arf_get_d(arb_midref(u), ARF_RND_NEAR) + 0.5);
    arb_sub_si(u, u, *kout, SDB_PREC);
    arb_get_mag(m, u);
    if (mag_cmp_2exp_si(m, -1) >= 0)
        goto cleanup;

    /* a = max(-Q, 0) + max(-s Re nu, 0), upper bounds */
    arb_neg(u, Q);
    arb_get_ubound_arf(arb_midref(a), u, SDB_PREC);
    if (arf_sgn(arb_midref(a)) < 0)
        arf_zero(arb_midref(a));
    if (s == 1)
        arb_neg(u, acb_realref(nu));
    else
        arb_set(u, acb_realref(nu));
    arb_get_ubound_arf(arb_midref(u), u, SDB_PREC);
    mag_zero(arb_radref(u));
    if (arf_sgn(arb_midref(u)) > 0)
        arb_add(a, a, u, SDB_PREC);
    arb_get_ubound_arf(arb_midref(a), a, SDB_PREC);
    mag_zero(arb_radref(a));
    if (!arb_is_finite(a))
        goto cleanup;

    /* e^G */
    _sd_bound_phi(t, e, e, NULL, alpha, sym, nu, w);
    arb_sub(G, acb_realref(t), acb_realref(c), SDB_PREC);
    sdb_exp_upper(eg, G);

    /* the weights at we */
    sdb_weights(g, m, sym, w, w);
    mag_add_ui(m, m, 1);

    if (nw > 0)
    {
        arb_sub_ui(u, a, 1, SDB_PREC);
        arb_get_ubound_arf(arb_midref(u), u, SDB_PREC);
        mag_zero(arb_radref(u));
        sdb_gamma_tail(res, K, u);
        mag_mul(res, res, eg);
    }
    if (nw > 1)
    {
        sdb_gamma_tail(res + 1, K, a);
        mag_mul(res + 1, res + 1, eg);
        if (nw > 2)
            mag_mul(res + 2, res + 1, m);
        mag_mul(res + 1, res + 1, g);
    }

    ok = 1;

cleanup:
    acb_clear(w); acb_clear(e); acb_clear(t);
    arb_clear(P); arb_clear(Q); arb_clear(K); arb_clear(a); arb_clear(G); arb_clear(u); arb_clear(v);
    mag_clear(eg); mag_clear(m); mag_clear(g);
    return ok;
}

/* The ray w = we + t d, t >= 0, to the linear valley.  If Re d <= 0
   (and Re d = 0 if beta != 0), then |alpha e^w| + |beta e^-w| <=
   |alpha| e^xe + |beta| e^-xe on the ray (xe = Re we), and if
   Re(nu d) >= mu > 0,

       Re(phi(w) - c) <= G - mu t,
       G = |alpha| e^xe + |beta| e^-xe - Re(nu we) - Re c,

   |dw| = |d| dt, g_1(w) <= g_1(we) and |w| <= |we| + |d| t. */
static int
sdb_tail_linear(mag_ptr res, int nw, const acb_t alpha, int sym, const acb_t nu,
    cd we, const acb_t c, cd d)
{
    acb_t w, dd, t;
    arb_t mu, G, u;
    mag_t eg, m, g, ad, imu;
    int ok = 0;

    if (!(d.re <= 0.0) || (sym && d.re != 0.0) || !(fabs(d.re) + fabs(d.im) < 1e300))
        return 0;

    acb_init(w); acb_init(dd); acb_init(t);
    arb_init(mu); arb_init(G); arb_init(u);
    mag_init(eg); mag_init(m); mag_init(g); mag_init(ad); mag_init(imu);

    arb_set_d(acb_realref(w), we.re);
    arb_set_d(acb_imagref(w), we.im);
    arb_set_d(acb_realref(dd), d.re);
    arb_set_d(acb_imagref(dd), d.im);

    /* mu = a lower bound for Re(nu d) */
    acb_mul(t, nu, dd, SDB_PREC);
    arb_get_lbound_arf(arb_midref(mu), acb_realref(t), SDB_PREC);
    mag_zero(arb_radref(mu));
    if (!arb_is_finite(mu) || !arb_is_positive(mu))
        goto cleanup;

    /* G (in arb arithmetic: the terms cancel) */
    sdb_weights(g, m, sym, w, w);
    arb_exp(G, acb_realref(w), SDB_PREC);
    if (sym)
    {
        arb_inv(u, G, SDB_PREC);
        arb_add(G, G, u, SDB_PREC);
    }
    acb_abs(u, alpha, SDB_PREC);
    arb_mul(G, G, u, SDB_PREC);             /* |alpha| (e^xe + [sym] e^-xe) */
    acb_mul(t, nu, w, SDB_PREC);
    arb_sub(G, G, acb_realref(t), SDB_PREC);
    arb_sub(G, G, acb_realref(c), SDB_PREC);
    sdb_exp_upper(eg, G);

    acb_get_mag(ad, dd);
    arb_inv(u, mu, SDB_PREC);
    arb_get_mag(imu, u);                    /* 1 / mu */

    /* |d| / mu, g_1 |d| / mu, |we| |d| / mu + |d|^2 / mu^2 */
    mag_mul(imu, imu, ad);
    mag_mul(eg, eg, imu);
    if (nw > 0)
        mag_set(res, eg);
    if (nw > 1)
        mag_mul(res + 1, eg, g);
    if (nw > 2)
    {
        mag_add(m, m, imu);
        mag_mul(res + 2, eg, m);
    }

    ok = 1;

cleanup:
    acb_clear(w); acb_clear(dd); acb_clear(t);
    arb_clear(mu); arb_clear(G); arb_clear(u);
    mag_clear(eg); mag_clear(m); mag_clear(g); mag_clear(ad); mag_clear(imu);
    return ok;
}

int
_sd_bound_tail(mag_ptr res, int nw, const acb_t alpha, int sym, const acb_t nu,
    cd we, const acb_t c, sd_end end, cd lin_dir)
{
    int k;

    if (!(fabs(we.re) + fabs(we.im) < 1e300))
        return 0;

    if (end.type == SD_PLUS || end.type == SD_MINUS)
        return sdb_tail_exp(res, nw, &k, alpha, sym, nu, we, c, (end.type == SD_PLUS) ? 1 : -1)
            && k == end.k;

    if (end.type == SD_LINEAR)
        return sdb_tail_linear(res, nw, alpha, sym, nu, we, c, lin_dir);

    return 0;
}

/* res[k] += bounds for the integrals over the polygon from the vertex
   i to the vertex j > i.  Any path between the two vertices will do:
   the straight segment is used if its bound (for the first weight) is
   below exp(logtol), and otherwise the range is split.  If P is not
   NULL, segments whose bounds are estimated in double precision (with
   the reference level top) to be too large are not tried. */
static void
sdb_polygon(mag_ptr res, int nw, const acb_t alpha, int sym, const acb_t nu,
    const cd * pts, const acb_t c, int i, int j, const mag_t tol,
    const sd_problem * P, double top, double logtol)
{
    mag_struct t[3];
    acb_t a, b;
    double g[2];
    int k, split;

    acb_init(a); acb_init(b);
    for (k = 0; k < 3; k++)
        mag_init(t + k);

    split = (j > i + 1) && (tol == NULL
        || (P != NULL && !(_sd_panel(g, P, pts[i], pts[j], top) <= logtol)));

    if (!split)
    {
        acb_set_d_d(a, pts[i].re, pts[i].im);
        acb_set_d_d(b, pts[j].re, pts[j].im);
        _sd_bound_segment(t, nw, alpha, sym, nu, a, b, c);
        split = (j > i + 1) && mag_cmp(t, tol) > 0;
    }

    if (!split)
    {
        for (k = 0; k < nw; k++)
            mag_add(res + k, res + k, t + k);
    }
    else
    {
        int m = i + (j - i) / 2;
        sdb_polygon(res, nw, alpha, sym, nu, pts, c, i, m, tol, P, top, logtol);
        sdb_polygon(res, nw, alpha, sym, nu, pts, c, m, j, tol, P, top, logtol);
    }

    acb_clear(a); acb_clear(b);
    for (k = 0; k < 3; k++)
        mag_clear(t + k);
}

/* The integral of |exp(phi - c)| |dw| from the vertex j of the piece to
   its end (the first if end = 0, the last otherwise) and along the ray
   from there into the valley, if the piece ends in a valley.  The
   polygon is replaced by fewer segments where their bounds stay below
   exp(logtol), which only serves to choose the segments, as does P (see
   above; it can be NULL).  Returns 0 if the ray is not certified. */
int
_sd_bound_end(mag_t res, const acb_t alpha, int sym, const acb_t nu,
    const sd_piece * pc, const acb_t c, int j, int end, cd lin_dir,
    const sd_problem * P, double logtol)
{
    int last = end ? pc->npts - 1 : 0, ok = 1;
    int valley = (pc->end[end].type == SD_PLUS || pc->end[end].type == SD_MINUS
        || pc->end[end].type == SD_LINEAR);

    mag_zero(res);
    if (valley)
        ok = _sd_bound_tail(res, 1, alpha, sym, nu, pc->pts[last], c, pc->end[end], lin_dir);

    if (ok && j != last)
    {
        mag_t tol;
        mag_init(tol);
        _sd_mag_approx_exp_d(tol, logtol);
        sdb_polygon(res, 1, alpha, sym, nu, pc->pts, c, FLINT_MIN(j, last), FLINT_MAX(j, last),
            tol, P, pc->c.re, logtol - 1.0);
        mag_clear(tol);
    }

    return ok && mag_is_finite(res);
}

/* The weighted integrals of |exp(phi - c)| |dw| between the vertices i
   and j > i of a polygon, added to res, with the polygon replaced by
   fewer straight segments where that does not cost much: the segment
   between the two vertices is used if its bound is not much larger
   than the double precision estimate of the bounds for the polygon
   between them (m[j] - m[i], where m[k] is the sum of the estimates for
   the first k edges), or small compared to the estimate for the whole
   polygon (m[n]); otherwise the range is split.  The estimates only
   serve to choose the segments. */
static void
sdb_polygon_merged(mag_ptr res, const acb_t alpha, int sym, const acb_t nu,
    const cd * pts, const acb_t c, int i, int j, const sd_problem * P,
    double top, const double * m, int n)
{
    mag_struct t[3];
    mag_t tol;
    acb_t a, b;
    double g[2], thr;
    int k, split;

    acb_init(a); acb_init(b);
    mag_init(tol);
    for (k = 0; k < 3; k++)
        mag_init(t + k);

    thr = 1.5 * (m[j] - m[i]) + 0.03 * m[n];

    split = (j > i + 1) && !(_sd_panel(g, P, pts[i], pts[j], top) <= log(thr));

    if (!split)
    {
        acb_set_d_d(a, pts[i].re, pts[i].im);
        acb_set_d_d(b, pts[j].re, pts[j].im);
        _sd_bound_segment(t, 3, alpha, sym, nu, a, b, c);
        if (j > i + 1)
        {
            mag_set_d(tol, 1.25 * thr);
            split = mag_cmp(t, tol) > 0;
        }
    }

    if (!split)
    {
        for (k = 0; k < 3; k++)
            mag_add(res + k, res + k, t + k);
    }
    else
    {
        int mid = i + (j - i) / 2;
        sdb_polygon_merged(res, alpha, sym, nu, pts, c, i, mid, P, top, m, n);
        sdb_polygon_merged(res, alpha, sym, nu, pts, c, mid, j, P, top, m, n);
    }

    acb_clear(a); acb_clear(b);
    mag_clear(tol);
    for (k = 0; k < 3; k++)
        mag_clear(t + k);
}

/* The weighted L1 norms of exp(phi) along the contour, as absolute
   bounds (not relative to the reference values of the pieces).  With a
   starting point (a ball), the contour begins with the segment from
   there to the first vertex of the first piece. */
int
_sd_contour_bounds(mag_ptr L, const sd_contour * C, const sd_problem * P,
    const acb_t alpha, int sym, const acb_t nu, const acb_t start)
{
    mag_struct S[3], G[3], t[3];
    acb_t a, b, c, q;
    arb_t u;
    mag_t f, f2;
    int i, k, j, ok = 1, any_geom = 0;

    if (C->npieces == 0)
        return 0;

    acb_init(a); acb_init(b); acb_init(c); acb_init(q);
    arb_init(u);
    mag_init(f); mag_init(f2);
    for (k = 0; k < 3; k++)
    {
        mag_init(S + k); mag_init(G + k); mag_init(t + k);
        mag_zero(L + k);
    }

    for (i = 0; i < C->npieces && ok; i++)
    {
        const sd_piece * pc = C->piece + i;

        acb_set_d_d(c, pc->c.re, pc->c.im);
        for (k = 0; k < 3; k++)
            mag_zero(S + k);

        if (i == 0 && start != NULL)
        {
            acb_set_d_d(b, pc->pts[0].re, pc->pts[0].im);
            _sd_bound_segment(t, 3, alpha, sym, nu, start, b, c);
            for (k = 0; k < 3; k++)
                mag_add(S + k, S + k, t + k);
        }

        /* the polygon, with fewer segments where the bounds allow it
           (any path between two vertices will do) */
        if (pc->npts >= 2)
        {
            double * m = flint_malloc(sizeof(double) * pc->npts);
            double g[2];

            m[0] = 0.0;
            for (j = 0; j + 1 < pc->npts; j++)
                m[j + 1] = m[j] + exp(_sd_panel(g, P, pc->pts[j], pc->pts[j + 1], pc->c.re));

            if (m[pc->npts - 1] > 0.0 && m[pc->npts - 1] < 1e300)
                sdb_polygon_merged(S, alpha, sym, nu, pc->pts, c, 0, pc->npts - 1,
                    P, pc->c.re, m, pc->npts - 1);
            else
            {
                for (j = 0; j + 1 < pc->npts; j++)
                {
                    acb_set_d_d(a, pc->pts[j].re, pc->pts[j].im);
                    acb_set_d_d(b, pc->pts[j + 1].re, pc->pts[j + 1].im);
                    _sd_bound_segment(t, 3, alpha, sym, nu, a, b, c);
                    for (k = 0; k < 3; k++)
                        mag_add(S + k, S + k, t + k);
                }
            }

            flint_free(m);
        }

        for (j = 0; j < 2 && ok; j++)
        {
            if (pc->end[j].type == SD_PLUS || pc->end[j].type == SD_MINUS
                || pc->end[j].type == SD_LINEAR)
            {
                ok = _sd_bound_tail(t, 3, alpha, sym, nu,
                    pc->pts[j ? pc->npts - 1 : 0], c, pc->end[j], P->lin_dir);
                for (k = 0; k < 3; k++)
                    mag_add(S + k, S + k, t + k);
            }
        }

        /* times exp(Re c) */
        arb_set_d(u, pc->c.re);
        sdb_exp_upper(f, u);
        for (k = 0; k < 3; k++)
        {
            mag_mul(S + k, S + k, f);
            if (pc->geom)
                mag_add(G + k, G + k, S + k);
            else
                mag_add(L + k, L + k, S + k);
        }
        any_geom |= pc->geom;
    }

    /* The pieces divided by 1 - q, q = exp(-2 pi i sigma nu): the
       derivative of I / (1 - q) with respect to nu is
       I' / (1 - q) + I q' / (1 - q)^2, |q'| = 2 pi |q|. */
    if (ok && any_geom)
    {
        acb_const_pi(q, SDB_PREC);
        acb_mul(q, q, nu, SDB_PREC);
        acb_mul_onei(q, q);
        acb_mul_si(q, q, (P->lin_dir.im > 0.0) ? -2 : 2, SDB_PREC);
        acb_exp(q, q, SDB_PREC);
        /* (the geometric series must converge: |q| < 1) */
        acb_get_mag(f, q);
        ok = mag_cmp_2exp_si(f, 0) < 0;
        acb_sub_ui(a, q, 1, SDB_PREC);
        acb_inv(a, a, SDB_PREC);
        acb_get_mag(f, a);                  /* 1 / |1 - q| */
        acb_get_mag(f2, q);
        mag_mul(f2, f2, f);
        mag_mul(f2, f2, f);
        mag_const_pi(t);
        mag_mul_2exp_si(t, t, 1);
        mag_mul(f2, f2, t);                 /* 2 pi |q| / |1 - q|^2 */

        mag_addmul(L + 2, G, f2);
        for (k = 0; k < 3; k++)
            mag_addmul(L + k, G + k, f);
    }

    for (k = 0; k < 3; k++)
    {
        if (!ok)
            mag_inf(L + k);
        mag_clear(S + k); mag_clear(G + k); mag_clear(t + k);
    }
    acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(q);
    arb_clear(u);
    mag_clear(f); mag_clear(f2);
    return ok;
}
