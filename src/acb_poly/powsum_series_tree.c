/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "thread_support.h"
#include "arb_poly.h"
#include "acb_poly.h"

/*
    Computes sum_{k=0}^{n-1} q^k (a+k)^(-(s+x)) mod x^len using a
    transposed multipoint evaluation.

    With L_k = log(a+k), w_k = q^k (a+k)^(-s) and an exact centre c,

        sum_k w_k exp(-L_k x) = exp(-c x) sum_j (p_j / j!) x^j,

    where p_j = sum_k w_k b_k^j are weighted power sums of the points
    b_k = c - L_k. The generating function of the p_j is the rational
    function

        sum_k w_k / (1 - b_k x) = P(x) / Q(x),   Q = prod_k (1 - b_k x),

    which we build with a subproduct tree (sum of fractions), truncating
    every node mod x^len, followed by a single power series division.
    This is the transpose of multipoint evaluation of a polynomial of
    length len at the points b_k. Over (n, len, prec) = Theta(m) the cost
    is O~(m^2) bit operations instead of O(m^3) for the naive sum.

    The tree is numerically unstable (the coefficients of Q can be
    exponentially large compared to the result), so we work with
    O(min(n, len)) guard bits, and we always work with exact (midpoint)
    inputs: radii of s, a, q are accounted for afterwards by a
    term-by-term perturbation bound, since feeding balls into the tree
    would amplify their radii in the same way as rounding errors.
*/

/* Power series with real part re and optional imaginary part im
   (im == NULL means identically zero). */

typedef struct
{
    arb_ptr Pr;
    arb_ptr Pi;
    arb_ptr Qr;
    arb_ptr Qi;
    slong Plen;
    slong Qlen;
    slong m;            /* number of leaves */
}
pq_struct;

typedef struct
{
    const acb_struct * s;   /* exact */
    const acb_struct * a;   /* exact, Re(a) > 0 */
    const acb_struct * q;   /* exact */
    acb_struct c;           /* exact centre */
    int a_real;             /* points (and Q) real */
    int w_real;             /* weights (and P) real */
    int q_one;
    int a_ui;               /* a is a small positive integer */
    ulong a_val;
    int s_si;               /* s is a small integer */
    slong s_val;
    slong len;
    slong wp;
    mag_ptr Lmag;           /* upper bounds for |log(a+k)|, or NULL */
    mag_ptr Wmag;           /* upper bounds for |w_k|, or NULL */
}
pq_args_struct;

static void
pq_init(pq_struct * x, pq_args_struct * args)
{
    x->Pr = x->Pi = x->Qr = x->Qi = NULL;
    x->Plen = x->Qlen = x->m = 0;
}

static void
pq_clear(pq_struct * x, pq_args_struct * args)
{
    if (x->Pr != NULL) _arb_vec_clear(x->Pr, x->Plen);
    if (x->Pi != NULL) _arb_vec_clear(x->Pi, x->Plen);
    if (x->Qr != NULL) _arb_vec_clear(x->Qr, x->Qlen);
    if (x->Qi != NULL) _arb_vec_clear(x->Qi, x->Qlen);
}

static void
pq_alloc(pq_struct * x, slong Plen, slong Qlen, const pq_args_struct * args)
{
    x->Plen = Plen;
    x->Qlen = Qlen;
    x->Pr = _arb_vec_init(Plen);
    x->Pi = args->w_real ? NULL : _arb_vec_init(Plen);
    x->Qr = _arb_vec_init(Qlen);
    x->Qi = args->a_real ? NULL : _arb_vec_init(Qlen);
}

/* {R, n} = {A, Alen} * {B, Blen} mod x^n, zero-padded; R is not aliased */
static void
_mullow(arb_ptr R, arb_srcptr A, slong Alen, arb_srcptr B, slong Blen, slong n, slong prec)
{
    slong m;

    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);
    m = FLINT_MIN(n, Alen + Blen - 1);

    if (Alen >= Blen)
        _arb_poly_mullow(R, A, Alen, B, Blen, m, prec);
    else
        _arb_poly_mullow(R, B, Blen, A, Alen, m, prec);

    _arb_vec_zero(R + m, n - m);
}

/* (Rr + i Ri) = (Ar + i Ai) (Br + i Bi) mod x^n; imaginary parts may be
   NULL (zero); Ri must be non-NULL unless both Ai and Bi are NULL.
   If add is set, the product is added to R instead. */
static void
_cmullow(arb_ptr Rr, arb_ptr Ri, arb_srcptr Ar, arb_srcptr Ai, slong Alen,
    arb_srcptr Br, arb_srcptr Bi, slong Blen, slong n, int add, slong prec)
{
    arb_ptr t;

    t = _arb_vec_init(n);

    _mullow(t, Ar, Alen, Br, Blen, n, prec);
    if (Ai != NULL && Bi != NULL)
    {
        arb_ptr u = _arb_vec_init(n);
        _mullow(u, Ai, Alen, Bi, Blen, n, prec);
        _arb_vec_sub(t, t, u, n, prec);
        _arb_vec_clear(u, n);
    }

    if (add)
        _arb_vec_add(Rr, Rr, t, n, prec);
    else
        _arb_vec_swap(Rr, t, n);

    if (Ai != NULL || Bi != NULL)
    {
        if (Ai == NULL)
            _mullow(t, Ar, Alen, Bi, Blen, n, prec);
        else if (Bi == NULL)
            _mullow(t, Ai, Alen, Br, Blen, n, prec);
        else
        {
            arb_ptr u = _arb_vec_init(n);
            _mullow(t, Ar, Alen, Bi, Blen, n, prec);
            _mullow(u, Ai, Alen, Br, Blen, n, prec);
            _arb_vec_add(t, t, u, n, prec);
            _arb_vec_clear(u, n);
        }

        if (add)
            _arb_vec_add(Ri, Ri, t, n, prec);
        else
            _arb_vec_swap(Ri, t, n);
    }

    _arb_vec_clear(t, n);
}

/* X <- X * (1 - b x) mod x^len, where X has length Xlen <= len and the
   result has length min(Xlen + 1, len). The caller provides space. */
static void
_mul_linear(arb_ptr Xr, arb_ptr Xi, slong Xlen, const arb_t br, const arb_t bi,
    slong len, slong prec)
{
    slong j;

    if (Xlen < len)
    {
        arb_zero(Xr + Xlen);
        if (Xi != NULL)
            arb_zero(Xi + Xlen);
        Xlen++;
    }

    for (j = Xlen - 1; j >= 1; j--)
    {
        if (bi == NULL)
        {
            arb_submul(Xr + j, br, Xr + j - 1, prec);
            if (Xi != NULL)
                arb_submul(Xi + j, br, Xi + j - 1, prec);
        }
        else
        {
            /* (Xr + i Xi)[j] -= (br + i bi) (Xr + i Xi)[j-1] */
            arb_submul(Xr + j, br, Xr + j - 1, prec);
            arb_addmul(Xr + j, bi, Xi + j - 1, prec);
            arb_submul(Xi + j, br, Xi + j - 1, prec);
            arb_submul(Xi + j, bi, Xr + j - 1, prec);
        }
    }
}

/* leaves lo <= k < hi, accumulated as P/Q += w_k / (1 - b_k x) */
static void
pq_basecase(pq_struct * res, slong lo, slong hi, pq_args_struct * args)
{
    slong k, j, len, wp, m, Plen, Qlen;
    acb_t ak, L, w, b, qpow, t;

    len = args->len;
    wp = args->wp;
    m = hi - lo;

    pq_alloc(res, FLINT_MIN(m, len), FLINT_MIN(m + 1, len), args);
    res->m = m;

    acb_init(ak);
    acb_init(L);
    acb_init(w);
    acb_init(b);
    acb_init(qpow);
    acb_init(t);

    if (!args->q_one)
        acb_pow_ui(qpow, args->q, lo, wp);

    Plen = Qlen = 0;

    for (k = lo; k < hi; k++)
    {
        /* L = log(a + k) */
        if (args->a_ui)
        {
            if (k == lo || args->a_val + k - 1 <= 1)
                arb_log_ui(acb_realref(L), args->a_val + k, wp);
            else
                arb_log_ui_from_prev(acb_realref(L), args->a_val + k,
                    acb_realref(L), args->a_val + k - 1, wp);
            arb_zero(acb_imagref(L));
        }
        else
        {
            acb_add_ui(ak, args->a, k, wp);
            acb_log(L, ak, wp);
        }

        /* w = q^k (a + k)^(-s) */
        if (args->s_si)
        {
            if (args->a_ui)
            {
                acb_set_ui(ak, args->a_val + k);
                acb_pow_si(w, ak, -args->s_val, wp);
            }
            else
            {
                acb_add_ui(ak, args->a, k, wp);
                acb_pow_si(w, ak, -args->s_val, wp);
            }
        }
        else
        {
            acb_mul(w, L, args->s, wp);
            acb_neg(w, w);
            acb_exp(w, w, wp);
        }

        if (!args->q_one)
        {
            acb_mul(w, w, qpow, wp);
            if (k < hi - 1)
                acb_mul(qpow, qpow, args->q, wp);
        }

        if (args->a_real)
            arb_zero(acb_imagref(L));
        if (args->w_real)
            arb_zero(acb_imagref(w));

        if (args->Lmag != NULL)
        {
            acb_get_mag(args->Lmag + k, L);
            acb_get_mag(args->Wmag + k, w);
        }

        /* b = c - L */
        acb_sub(b, &args->c, L, wp);

        if (k == lo)
        {
            arb_set(res->Pr, acb_realref(w));
            if (res->Pi != NULL)
                arb_set(res->Pi, acb_imagref(w));
            Plen = 1;

            arb_one(res->Qr);
            if (res->Qi != NULL)
                arb_zero(res->Qi);
            Qlen = 1;
        }
        else
        {
            /* P <- P (1 - b x) + w Q */
            _mul_linear(res->Pr, res->Pi, Plen, acb_realref(b),
                args->a_real ? NULL : acb_imagref(b), len, wp);
            Plen = FLINT_MIN(Plen + 1, len);

            for (j = 0; j < Plen; j++)
            {
                /* Q has length >= Plen here */
                if (res->Qi == NULL)
                {
                    arb_addmul(res->Pr + j, acb_realref(w), res->Qr + j, wp);
                    if (res->Pi != NULL)
                        arb_addmul(res->Pi + j, acb_imagref(w), res->Qr + j, wp);
                }
                else
                {
                    acb_set_arb_arb(t, res->Qr + j, res->Qi + j);
                    acb_mul(t, t, w, wp);
                    arb_add(res->Pr + j, res->Pr + j, acb_realref(t), wp);
                    arb_add(res->Pi + j, res->Pi + j, acb_imagref(t), wp);
                }
            }
        }

        /* Q <- Q (1 - b x) */
        _mul_linear(res->Qr, res->Qi, Qlen, acb_realref(b),
            args->a_real ? NULL : acb_imagref(b), len, wp);
        Qlen = FLINT_MIN(Qlen + 1, len);
    }

    FLINT_ASSERT(Plen == res->Plen);
    FLINT_ASSERT(Qlen == res->Qlen);

    acb_clear(ak);
    acb_clear(L);
    acb_clear(w);
    acb_clear(b);
    acb_clear(qpow);
    acb_clear(t);
}

static void
pq_merge(pq_struct * res, pq_struct * left, pq_struct * right, pq_args_struct * args)
{
    slong len, wp, Plen, Qlen;

    len = args->len;
    wp = args->wp;

    /* a node with m leaves has deg P <= m - 1, deg Q <= m */
    res->m = left->m + right->m;
    Plen = FLINT_MIN(len, res->m);
    Qlen = FLINT_MIN(len, res->m + 1);

    pq_alloc(res, Plen, Qlen, args);

    /* P = P1 Q2 + P2 Q1 */
    _cmullow(res->Pr, res->Pi, left->Pr, left->Pi, left->Plen,
        right->Qr, right->Qi, right->Qlen, Plen, 0, wp);
    _cmullow(res->Pr, res->Pi, right->Pr, right->Pi, right->Plen,
        left->Qr, left->Qi, left->Qlen, Plen, 1, wp);

    /* Q = Q1 Q2 */
    _cmullow(res->Qr, res->Qi, left->Qr, left->Qi, left->Qlen,
        right->Qr, right->Qi, right->Qlen, Qlen, 0, wp);
}

/* Upper bound perturbation term for inexact inputs; see the comments
   in the main function. */
static void
_powsum_perturbation(mag_ptr E, const acb_t s, const acb_t a, const acb_t q,
    const acb_t sm, const acb_t am, const acb_t qm,
    mag_srcptr Lmag, mag_srcptr Wmag, slong n, slong len)
{
    mag_t rs, ra, rq, abs_s, dq, dl, lam, eta, t, u, lammax, lamfloor;
    mag_ptr Lam, alpha, beta, BL, BA, BB;
    slong k, i, j, nb;
    acb_t ak;
    double logfloor, logstep;

    mag_init(rs); mag_init(ra); mag_init(rq); mag_init(abs_s);
    mag_init(dq); mag_init(dl); mag_init(lam); mag_init(eta);
    mag_init(t); mag_init(u); mag_init(lammax); mag_init(lamfloor);
    acb_init(ak);

    Lam = _mag_vec_init(n);
    alpha = _mag_vec_init(n);
    beta = _mag_vec_init(n);

    mag_add(rs, arb_radref(acb_realref(s)), arb_radref(acb_imagref(s)));
    mag_add(ra, arb_radref(acb_realref(a)), arb_radref(acb_imagref(a)));
    mag_add(rq, arb_radref(acb_realref(q)), arb_radref(acb_imagref(q)));
    acb_get_mag(abs_s, sm);

    /* dq = r_q / (|q_m| - r_q) */
    if (mag_is_zero(rq))
    {
        mag_zero(dq);
    }
    else
    {
        acb_get_mag_lower(t, qm);
        mag_sub_lower(t, t, rq);
        mag_div(dq, rq, t);
    }

    mag_zero(lammax);

    for (k = 0; k < n; k++)
    {
        /* dl = r_a / (|a_m + k| - r_a) bounds |log(a+k) - log(a_m+k)| */
        if (mag_is_zero(ra))
        {
            mag_zero(dl);
        }
        else
        {
            acb_add_ui(ak, am, k, MAG_BITS);
            acb_get_mag_lower(t, ak);
            mag_sub_lower(t, t, ra);
            mag_div(dl, ra, t);
        }

        /* Lam >= |log(a+k)| on the ball */
        mag_add(Lam + k, Lmag + k, dl);

        /* eta >= |log(u) - log(u_m)| with u = q^k (a+k)^(-s) */
        mag_mul_ui(eta, dq, k);
        mag_addmul(eta, rs, Lam + k);
        mag_addmul(eta, abs_s, dl);

        /* alpha = |u - u_m| <= |u_m| (exp(eta) - 1), beta = |u_m| dl */
        mag_expm1(t, eta);
        mag_mul(alpha + k, Wmag + k, t);
        mag_mul(beta + k, Wmag + k, dl);

        if (mag_cmp(Lam + k, lammax) > 0)
            mag_set(lammax, Lam + k);
    }

    /*
        E_i = (1/i!) sum_k (alpha_k Lam_k^i + i beta_k Lam_k^(i-1)).

        Doing this term by term costs n len operations, so we group the
        points into nb buckets (geometrically spaced between lammax/16
        and lammax) and use the largest Lam_k in each bucket. Bucket
        assignment is a heuristic; the representative is a rigorous max.
        The overestimation is at most a factor 16^(i/nb).
    */
    nb = FLINT_MIN(n, FLINT_MAX(len, 32));

    BL = _mag_vec_init(nb);
    BA = _mag_vec_init(nb);
    BB = _mag_vec_init(nb);

    mag_mul_2exp_si(lamfloor, lammax, -4);
    logfloor = mag_get_d_log2_approx(lamfloor);
    logstep = 4.0 / nb;

    for (k = 0; k < n; k++)
    {
        if (nb == n)
            j = k;
        else if (mag_cmp(Lam + k, lamfloor) <= 0)
            j = 0;
        else
        {
            double v = (mag_get_d_log2_approx(Lam + k) - logfloor) / logstep;
            j = (v <= 0.0) ? 0 : (slong) v;
            j = FLINT_MIN(j, nb - 1);
            j = FLINT_MAX(j, 0);
        }

        mag_max(BL + j, BL + j, Lam + k);
        mag_add(BA + j, BA + j, alpha + k);
        mag_add(BB + j, BB + j, beta + k);
    }

    for (i = 0; i < len; i++)
        mag_zero(E + i);

    for (j = 0; j < nb; j++)
    {
        if (mag_is_zero(BA + j) && mag_is_zero(BB + j))
            continue;

        /* t = BA * BL^i, u = BB * BL^(i-1) */
        mag_set(t, BA + j);
        mag_set(u, BB + j);

        mag_add(E, E, t);

        for (i = 1; i < len; i++)
        {
            mag_mul(t, t, BL + j);
            mag_add(E + i, E + i, t);
            mag_mul_ui(lam, u, i);
            mag_add(E + i, E + i, lam);
            mag_mul(u, u, BL + j);
        }
    }

    for (i = 2; i < len; i++)
    {
        mag_rfac_ui(t, i);
        mag_mul(E + i, E + i, t);
    }

    _mag_vec_clear(Lam, n);
    _mag_vec_clear(alpha, n);
    _mag_vec_clear(beta, n);
    _mag_vec_clear(BL, nb);
    _mag_vec_clear(BA, nb);
    _mag_vec_clear(BB, nb);

    mag_clear(rs); mag_clear(ra); mag_clear(rq); mag_clear(abs_s);
    mag_clear(dq); mag_clear(dl); mag_clear(lam); mag_clear(eta);
    mag_clear(t); mag_clear(u); mag_clear(lammax); mag_clear(lamfloor);
    acb_clear(ak);
}

/* core with explicit working precision; inputs s, a, q, c must be exact */
static void
_powsum_series_tree_exact(acb_ptr z, const acb_t s, const acb_t a,
    const acb_t q, const acb_t c, slong n, slong len, slong wp, slong prec,
    mag_ptr Lmag, mag_ptr Wmag)
{
    pq_args_struct args;
    pq_struct res;
    arb_ptr Qinv_r, Qinv_i, Rr, Ri, Er, Ei, Tr, Ti;
    slong i;

    args.s = s;
    args.a = a;
    args.q = q;
    args.a_real = acb_is_real(a);
    args.q_one = acb_is_one(q);
    args.w_real = args.a_real && acb_is_real(s) && acb_is_real(q);
    args.a_ui = args.a_real && arf_is_int(arb_midref(acb_realref(a)))
        && arf_cmp_2exp_si(arb_midref(acb_realref(a)), FLINT_BITS - 2) < 0;
    args.a_val = args.a_ui ? arf_get_si(arb_midref(acb_realref(a)), ARF_RND_DOWN) : 0;
    if (args.a_ui && (args.a_val < 1 || (ulong) n >= UWORD(1) << (FLINT_BITS - 2)))
        args.a_ui = 0;
    args.s_si = acb_is_real(s) && arf_is_int(arb_midref(acb_realref(s)))
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 20) < 0;
    args.s_val = args.s_si ? arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN) : 0;
    args.len = len;
    args.wp = wp;
    args.Lmag = Lmag;
    args.Wmag = Wmag;

    acb_init(&args.c);
    acb_set(&args.c, c);
    if (args.a_real)
        arb_zero(acb_imagref(&args.c));

    pq_init(&res, &args);

    flint_parallel_binary_splitting(&res,
        (bsplit_basecase_func_t) pq_basecase,
        (bsplit_merge_func_t) pq_merge,
        sizeof(pq_struct),
        (bsplit_init_func_t) pq_init,
        (bsplit_clear_func_t) pq_clear,
        &args, 0, n, 16, -1, 0);

    /* R = P / Q mod x^len */
    Qinv_r = _arb_vec_init(len);
    Qinv_i = args.a_real ? NULL : _arb_vec_init(len);
    Rr = _arb_vec_init(len);
    Ri = (args.a_real && args.w_real) ? NULL : _arb_vec_init(len);

    if (args.a_real)
    {
        _arb_poly_inv_series(Qinv_r, res.Qr, res.Qlen, len, wp);
    }
    else
    {
        acb_ptr Q, Qinv;
        Q = _acb_vec_init(res.Qlen);
        Qinv = _acb_vec_init(len);
        for (i = 0; i < res.Qlen; i++)
            acb_set_arb_arb(Q + i, res.Qr + i, res.Qi + i);
        _acb_poly_inv_series(Qinv, Q, res.Qlen, len, wp);
        for (i = 0; i < len; i++)
        {
            arb_swap(Qinv_r + i, acb_realref(Qinv + i));
            arb_swap(Qinv_i + i, acb_imagref(Qinv + i));
        }
        _acb_vec_clear(Q, res.Qlen);
        _acb_vec_clear(Qinv, len);
    }

    _cmullow(Rr, Ri, res.Pr, res.Pi, res.Plen, Qinv_r, Qinv_i, len, len, 0, wp);

    /* Borel transform: p_j -> p_j / j! */
    _arb_poly_borel_transform(Rr, Rr, len, wp);
    if (Ri != NULL)
        _arb_poly_borel_transform(Ri, Ri, len, wp);

    /* multiply by exp(-c x) */
    Er = _arb_vec_init(len);
    Ei = args.a_real ? NULL : _arb_vec_init(len);
    Tr = _arb_vec_init(len);
    Ti = (Ri == NULL) ? NULL : _arb_vec_init(len);

    arb_one(Er);
    if (Ei != NULL)
        arb_zero(Ei);
    for (i = 1; i < len; i++)
    {
        if (Ei == NULL)
        {
            arb_mul(Er + i, Er + i - 1, acb_realref(&args.c), wp);
        }
        else
        {
            acb_t t;
            acb_init(t);
            acb_set_arb_arb(t, Er + i - 1, Ei + i - 1);
            acb_mul(t, t, &args.c, wp);
            arb_swap(Er + i, acb_realref(t));
            arb_swap(Ei + i, acb_imagref(t));
            acb_clear(t);
        }
    }
    for (i = 1; i < len; i += 2)
    {
        arb_neg(Er + i, Er + i);
        if (Ei != NULL)
            arb_neg(Ei + i, Ei + i);
    }
    _arb_poly_borel_transform(Er, Er, len, wp);
    if (Ei != NULL)
        _arb_poly_borel_transform(Ei, Ei, len, wp);

    _cmullow(Tr, Ti, Rr, Ri, len, Er, Ei, len, len, 0, wp);

    for (i = 0; i < len; i++)
    {
        arb_set_round(acb_realref(z + i), Tr + i, prec);
        if (Ti != NULL)
            arb_set_round(acb_imagref(z + i), Ti + i, prec);
        else
            arb_zero(acb_imagref(z + i));
    }

    _arb_vec_clear(Qinv_r, len);
    if (Qinv_i != NULL) _arb_vec_clear(Qinv_i, len);
    _arb_vec_clear(Rr, len);
    if (Ri != NULL) _arb_vec_clear(Ri, len);
    _arb_vec_clear(Er, len);
    if (Ei != NULL) _arb_vec_clear(Ei, len);
    _arb_vec_clear(Tr, len);
    if (Ti != NULL) _arb_vec_clear(Ti, len);

    pq_clear(&res, &args);
    acb_clear(&args.c);
}

/* Exact centre c ~= mean of log(a+k), 0 <= k < n. Centring at the mean
   of the points (rather than the midpoint of their range) makes the
   linear coefficient of Q vanish and reduces the precision loss in the
   tree several times over. */
static void
_powsum_tree_centre(acb_t c, const acb_t a, slong n)
{
    acb_t t;
    acb_init(t);

    acb_add_ui(t, a, n, 64);
    acb_lgamma(t, t, 64);
    acb_lgamma(c, a, 64);
    acb_sub(c, t, c, 64);
    acb_div_ui(c, c, n, 64);

    if (!acb_is_finite(c))
    {
        acb_log(c, a, 64);
        acb_add_ui(t, a, n - 1, 64);
        acb_log(t, t, 64);
        acb_add(c, c, t, 64);
        acb_mul_2exp_si(c, c, -1);
    }

    acb_get_mid(c, c);
    arf_set_round(arb_midref(acb_realref(c)), arb_midref(acb_realref(c)), 32, ARF_RND_NEAR);
    arf_set_round(arb_midref(acb_imagref(c)), arb_midref(acb_imagref(c)), 32, ARF_RND_NEAR);

    if (!acb_is_finite(c))
        acb_zero(c);

    acb_clear(t);
}

/*
    Heuristic number of guard bits for the tree algorithm so that the
    output radii are comparable to those of the naive algorithm.

    The dominant loss comes from the final multiplication by exp(-c x):
    coefficient i is computed with absolute error proportional to
    sum_k |w_k| (|b_k| + |c|)^i / i! while its natural size (that
    obtained by the naive term-by-term summation) is
    sum_k |w_k| |log(a+k)|^i / i!. We estimate the maximum ratio in
    double precision over a sample of i and of the points k (which also
    captures the dependence on s, a, q through the weights), and add
    an empirical model for the remaining loss in the tree and the series
    division (fitted to measurements with n, len up to a few thousand).
*/
slong
_acb_poly_powsum_series_tree_guard_bits(const acb_t s, const acb_t a,
    const acb_t q, const acb_t c, slong n, slong len)
{
    double sr, si, ar, ai, lq, cr, ci, absc, proxy, resid, m;
    double *Lw, *LB, *LL;
    slong i, k, ns, S;

    m = FLINT_MIN(n, len);

    if (n >= len)
        resid = m * (0.7 + 0.45 * log2((double) n / m));
    else
        resid = m * (0.75 + 0.3 * log2((double) len / m));

    sr = arf_get_d(arb_midref(acb_realref(s)), ARF_RND_NEAR);
    si = arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_NEAR);
    ar = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    ai = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    cr = arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR);
    ci = arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR);
    lq = log(hypot(arf_get_d(arb_midref(acb_realref(q)), ARF_RND_NEAR),
                   arf_get_d(arb_midref(acb_imagref(q)), ARF_RND_NEAR)));
    absc = hypot(cr, ci);

    if (!(fabs(sr) < 1e15 && fabs(si) < 1e15 && fabs(ar) < 1e15 && fabs(ai) < 1e15
        && fabs(lq) < 1e10 && ar > 0.0))
    {
        /* should not happen in practice; something crude but safe-ish */
        return 2 * (n + len) + 64;
    }

    /* sample points: all k when n is small, otherwise the first few k
       individually and the rest in geometric strata with multiplicity */
    S = 960;
    ns = (n <= 1024) ? n : 64 + S;
    Lw = flint_malloc(sizeof(double) * ns);
    LB = flint_malloc(sizeof(double) * ns);
    LL = flint_malloc(sizeof(double) * ns);

    {
        slong klo, khi, idx = 0;
        double rho = 1.0;

        if (n > 1024)
            rho = pow((n - 1 + ar) / (64 + ar), 1.0 / S);

        klo = 0;
        while (klo < n)
        {
            double kr, lr, li;

            if (n <= 1024 || klo < 64)
                khi = klo + 1;
            else
            {
                khi = (slong) ((klo + ar) * rho - ar);
                khi = FLINT_MAX(khi, klo + 1);
                khi = FLINT_MIN(khi, n);
            }

            kr = 0.5 * (klo + khi - 1);
            lr = 0.5 * log((ar + kr) * (ar + kr) + ai * ai);
            li = atan2(ai, ar + kr);

            Lw[idx] = -(sr * lr - si * li) + kr * lq + log((double) (khi - klo));
            LB[idx] = log(hypot(cr - lr, ci - li) + absc);
            LL[idx] = (lr == 0.0 && li == 0.0) ? -1e300 : log(hypot(lr, li));
            idx++;
            klo = khi;

            if (idx >= ns)
                break;
        }

        ns = idx;
    }

    /* The loss for coefficient i is estimated as

            log2( (sum_k |w_k|) max_k (|b_k| + |c|)^i / sum_k |w_k| |log(a+k)|^i ).

       (The error in the tree is not localised to individual points, so
       we pair the total weight with the largest point.) Since
       |b_k| + |c| >= |log(a+k)|, this is nondecreasing in i, so it
       suffices to evaluate it at i = len - 1. */
    {
        double lsw, lmaxb, lr, li, B, t, mx;

        lmaxb = -1e300;
        for (k = 0; k < ns; k++)
            lmaxb = FLINT_MAX(lmaxb, LB[k]);

        /* make sure that the endpoints are included */
        lr = 0.5 * log((ar + n - 1) * (ar + n - 1) + ai * ai);
        li = atan2(ai, ar + n - 1);
        lmaxb = FLINT_MAX(lmaxb, log(hypot(cr - lr, ci - li) + absc));

        i = len - 1;

        /* log sum exp of Lw[k] */
        mx = -1e300;
        for (k = 0; k < ns; k++)
            mx = FLINT_MAX(mx, Lw[k]);
        t = 0.0;
        for (k = 0; k < ns; k++)
            t += exp(Lw[k] - mx);
        lsw = mx + log(t);

        /* log sum exp of Lw[k] + i LL[k] */
        mx = -1e300;
        for (k = 0; k < ns; k++)
            if (i == 0 || LL[k] > -1e300)
                mx = FLINT_MAX(mx, Lw[k] + i * LL[k]);
        t = 0.0;
        for (k = 0; k < ns; k++)
            if (i == 0 || LL[k] > -1e300)
                t += exp(Lw[k] + i * LL[k] - mx);
        B = mx + log(t);

        proxy = (lsw + i * lmaxb - B) / log(2.0);
        proxy = FLINT_MAX(proxy, 0.0);

        if (!(proxy < 1e15))
            proxy = 2 * (n + len);
    }

    flint_free(Lw);
    flint_free(LB);
    flint_free(LL);

    return (slong) (proxy + resid + 16.0);
}

void
_acb_poly_powsum_series_tree(acb_ptr z, const acb_t s, const acb_t a,
    const acb_t q, slong n, slong len, slong prec)
{
    acb_t sm, am, qm, c;
    mag_ptr Lmag, Wmag, E;
    slong i, wp;
    int exact, real;

    if (len <= 0)
        return;

    if (n <= 0)
    {
        _acb_vec_zero(z, len);
        return;
    }

    /* The perturbation bounds need Re(a) > 0 (no branch cut); q must
       be bounded away from zero when inexact. */
    if (!acb_is_finite(s) || !acb_is_finite(a) || !acb_is_finite(q) ||
        !arb_is_positive(acb_realref(a)) ||
        (!acb_is_exact(q) && acb_contains_zero(q)) || acb_is_zero(q))
    {
        _acb_poly_powsum_series_naive(z, s, a, q, n, len, prec);
        return;
    }

    acb_init(sm);
    acb_init(am);
    acb_init(qm);
    acb_init(c);

    acb_get_mid(sm, s);
    acb_get_mid(am, a);
    acb_get_mid(qm, q);

    exact = acb_is_exact(s) && acb_is_exact(a) && acb_is_exact(q);
    real = acb_is_real(s) && acb_is_real(a) && acb_is_real(q);

    _powsum_tree_centre(c, am, n);
    wp = prec + _acb_poly_powsum_series_tree_guard_bits(sm, am, qm, c, n, len);

    Lmag = Wmag = NULL;
    if (!exact)
    {
        Lmag = _mag_vec_init(n);
        Wmag = _mag_vec_init(n);
    }

    _powsum_series_tree_exact(z, sm, am, qm, c, n, len, wp, prec, Lmag, Wmag);

    if (!exact)
    {
        /*
            Let u_k = q^k (a+k)^(-s) and l_k = log(a+k), so that the
            coefficient of x^i is sum_k u_k (-l_k)^i / i!. For (s, a, q)
            in the input balls, with midpoints (s_m, a_m, q_m),

                |u l^i - u_m l_m^i| <= |u - u_m| Lam^i + |u_m| i |l - l_m| Lam^(i-1)

            where Lam bounds |l| on the ball. We bound
            |l - l_m| <= r_a / (|a_m + k| - r_a) and
            |log(u / u_m)| <= eta = k r_q / (|q_m| - r_q) + r_s Lam + |s_m| |l - l_m|,
            so that |u - u_m| <= |u_m| (exp(eta) - 1).
        */
        E = _mag_vec_init(len);
        _powsum_perturbation(E, s, a, q, sm, am, qm, Lmag, Wmag, n, len);

        for (i = 0; i < len; i++)
        {
            arb_add_error_mag(acb_realref(z + i), E + i);
            if (!real)
                arb_add_error_mag(acb_imagref(z + i), E + i);
        }

        _mag_vec_clear(E, len);
        _mag_vec_clear(Lmag, n);
        _mag_vec_clear(Wmag, n);
    }

    acb_clear(sm);
    acb_clear(am);
    acb_clear(qm);
    acb_clear(c);
}

/*
    Decide whether _acb_poly_powsum_series_tree is likely to be faster
    than _acb_poly_powsum_series_naive (or the threaded version, if
    naive_threaded is set). We compare simple cost models (in
    microseconds) fitted to timings on x86-64 for n, len up to a few
    thousand and prec up to 4096: the naive algorithm costs
    O(n len M(prec)) while the tree costs
    O((n log n + len log len) M(wp) / wp) with wp = prec + guard bits;
    the tree is only partially parallelised.
*/
int
_acb_poly_powsum_series_tree_is_faster(const acb_t s, const acb_t a,
    const acb_t q, slong n, slong len, slong prec, int naive_threaded)
{
    double N, L, w, cpx, tnaive, ttree, margin;
    slong wp, threads;
    acb_t sm, am, qm, c;

    if (n < 32 || len < 32 || n > WORD_MAX / 4 || len > WORD_MAX / 4)
        return 0;

    if (!acb_is_finite(s) || !acb_is_finite(a) || !acb_is_finite(q) ||
        !arb_is_positive(acb_realref(a)) || acb_contains_zero(q))
        return 0;

    N = n;
    L = len;
    cpx = !(acb_is_real(s) && acb_is_real(a) && acb_is_real(q));
    threads = flint_get_num_threads();

    w = prec / 64.0;
    tnaive = N * L * (0.183 + 0.00739 * pow(w, 1.458)) * (1.0 + 0.93 * cpx);
    if (naive_threaded)
        tnaive /= threads;

#define TREE_COST(w) ((0.231 * N * log2(2 * N) * pow(w, 0.957) \
                      + 0.140 * L * log2(2 * L) * pow(w, 0.957) \
                      + 0.0025 * N * pow(w, 2.317)) * (1.0 + 0.83 * cpx) \
                      * (0.6 + 0.4 / threads))

    /* the model was fitted with q = 1; with geometrically varying weights
       (e.g. polylogarithms) the tree is somewhat slower near the
       crossover, so require a larger predicted gain */
    margin = acb_is_one(q) ? 1.25 : 1.5;

    /* lower bound without guard bits */
    ttree = TREE_COST(w);
    if (tnaive < margin * ttree)
        return 0;

    acb_init(sm);
    acb_init(am);
    acb_init(qm);
    acb_init(c);

    acb_get_mid(sm, s);
    acb_get_mid(am, a);
    acb_get_mid(qm, q);

    _powsum_tree_centre(c, am, n);
    wp = prec + _acb_poly_powsum_series_tree_guard_bits(sm, am, qm, c, n, len);

    acb_clear(sm);
    acb_clear(am);
    acb_clear(qm);
    acb_clear(c);

    if (wp > 100 * prec + 100000)
        return 0;

    ttree = TREE_COST(wp / 64.0);

#undef TREE_COST

    return tnaive > margin * ttree;
}
