/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb_poly.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

static void
bsplit(acb_t A, acb_t B, acb_t C, acb_t D,
    const acb_t b, const acb_t z, slong n0, slong n1, slong prec)
{
    if (n1 - n0 == 1)
    {
        acb_zero(A);
        acb_one(B);
        acb_neg(C, b);
        acb_add_si(C, C, 2 - n0, prec);
        acb_mul_si(C, C, n0 - 1, prec);
        acb_sub(D, z, b, prec);
        acb_add_si(D, D, 2 - 2 * n0, prec);
    }
    else
    {
        slong m;
        acb_t T, A2, B2, C2, D2;

        acb_init(T);
        acb_init(A2);
        acb_init(B2);
        acb_init(C2);
        acb_init(D2);

        m = n0 + (n1 - n0) / 2;

        bsplit(A, B, C, D, b, z, n0, m, prec);
        bsplit(A2, B2, C2, D2, b, z, m, n1, prec);

        acb_set(T, A);
        acb_mul(A, A, A2, prec);
        acb_addmul(A, B2, C, prec);
        acb_mul(C, C, D2, prec);
        acb_addmul(C, C2, T, prec);

        acb_set(T, B);
        acb_mul(B, B, A2, prec);
        acb_addmul(B, B2, D, prec);
        acb_mul(D, D, D2, prec);
        acb_addmul(D, C2, T, prec);

        acb_clear(T);
        acb_clear(A2);
        acb_clear(B2);
        acb_clear(C2);
        acb_clear(D2);
    }
}

static void
acb_hypgeom_u_si_rec(acb_t res, slong a, const acb_t b, const acb_t z, slong prec)
{
    slong k;
    acb_t u0, u1, t;

    if (a > 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (a == 0)
    {
        acb_one(res);
        return;
    }
    else if (a == -1)
    {
        acb_sub(res, z, b, prec);
        return;
    }

    /* special-case U(-n, -n+1, z) = z^n */
    if (acb_equal_si(b, a + 1))
    {
        acb_pow_si(res, z, -a, prec);
        return;
    }

    acb_init(u0);
    acb_init(u1);
    acb_init(t);

    acb_one(u0);
    acb_sub(u1, z, b, prec);

    if (-a < 20)
    {
        for (k = 2; k <= -a; k++)
        {
            acb_neg(t, b);
            acb_add_si(t, t, 2 - k, prec);
            acb_mul_si(t, t, k - 1, prec);
            acb_mul(u0, u0, t, prec);

            acb_sub(t, z, b, prec);
            acb_add_si(t, t, 2 - 2 * k, prec);
            acb_addmul(u0, u1, t, prec);

            acb_swap(u0, u1);
        }

        acb_set(res, u1);
    }
    else
    {
        acb_t A, B, C, D;

        acb_init(A);
        acb_init(B);
        acb_init(C);
        acb_init(D);

        bsplit(A, B, C, D, b, z, 2, -a + 1, prec);

        acb_sub(A, z, b, prec);
        acb_mul(D, D, A, prec);
        acb_add(res, C, D, prec);

        acb_clear(A);
        acb_clear(B);
        acb_clear(C);
        acb_clear(D);
    }

    acb_clear(u0);
    acb_clear(u1);
    acb_clear(t);
}

void
acb_hypgeom_u_1f1_series(acb_poly_t res,
    const acb_poly_t a, const acb_poly_t b, const acb_poly_t z,
    slong len, slong prec)
{
    acb_poly_t s, u, A, B;
    acb_poly_struct aa[3];
    arb_t c;
    slong wlen;
    int singular;

    acb_poly_init(s);
    acb_poly_init(u);
    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(aa + 0);
    acb_poly_init(aa + 1);
    acb_poly_init(aa + 2);
    arb_init(c);

    singular = (b->length == 0) || acb_is_int(b->coeffs);
    wlen = len + (singular != 0);

    /* A = rgamma(a-b+1) * 1F~1(a,b,z) */
    acb_poly_sub(u, a, b, prec);
    acb_poly_add_si(u, u, 1, prec);
    acb_poly_rgamma_series(A, u, wlen, prec);

    /* todo: handle a = 1 efficiently */
    acb_poly_set(aa, a);
    acb_poly_set(aa + 1, b);
    acb_poly_one(aa + 2);
    acb_hypgeom_pfq_series_direct(s, aa, 1, aa + 1, 2, z, 1, -1, wlen, prec);
    acb_poly_mullow(A, A, s, wlen, prec);

    /* B = rgamma(a) * 1F~1(a-b+1,2-b,z) * z^(1-b) */
    acb_poly_set(aa, u);
    acb_poly_add_si(aa + 1, b, -2, prec);
    acb_poly_neg(aa + 1, aa + 1);
    acb_hypgeom_pfq_series_direct(s, aa, 1, aa + 1, 2, z, 1, -1, wlen, prec);
    acb_poly_rgamma_series(B, a, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_add_si(u, b, -1, prec);
    acb_poly_neg(u, u);
    acb_poly_pow_series(s, z, u, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_sub(A, A, B, prec);

    /* multiply by pi csc(pi b) */
    acb_poly_sin_pi_series(B, b, wlen, prec);

    if (singular)
    {
        acb_poly_shift_right(A, A, 1);
        acb_poly_shift_right(B, B, 1);
    }

    acb_poly_div_series(res, A, B, len, prec);

    arb_const_pi(c, prec);
    _acb_vec_scalar_mul_arb(res->coeffs, res->coeffs, res->length, c, prec);

    acb_poly_clear(s);
    acb_poly_clear(u);
    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(aa + 0);
    acb_poly_clear(aa + 1);
    acb_poly_clear(aa + 2);
    arb_clear(c);
}

/* ------------------------------------------------------------------------- */
/* b close to an integer                                                      */
/* ------------------------------------------------------------------------- */

/* For b = n + eps with a small eps, the two terms of the connection
   formula are of the size 1 / eps and cancel; on a ball for b of
   radius rad the result has log2(|eps| / rad) bits at best, and a ball
   which contains n gives nothing.  U is analytic in b:

       U(a, n + d, z) = (-1)^n pi (A(d) - B(d)) / sin(pi d),

       A = 1F~1(a, b, z) / Gamma(a - b + 1),
       B = z^(1-b) 1F~1(a - b + 1, 2 - b, z) / Gamma(a),

   where A and B are entire functions of d which agree at d = 0.  The
   Taylor coefficients u_j of U at d = 0 are those of the power series
   U(a, n + x, z) (acb_hypgeom_u_1f1_series), and by Cauchy's estimate
   on the circle |d| = R = 1/4, where |sin(pi d)| >= sin(pi/4),

       |u_j| <= M R^(-j),    M = pi sqrt(2) sup (|A| + |B|),

   so that the remainder after J terms is at most M q^J / (1 - q) with
   q = |eps| / R.  The supremum is bounded from the terms of the two
   series: with be in the ball with the centre be0 (an integer) and the
   radius R in both parts, and K = max(0, 2 - be0),

       |1F~1(al, be, z)| <= sum_(k<K) c_k g_k + c_K G T,
       c_k = sup |(al)_k| |z|^k / k!,
       G = sup |1 / Gamma(be + K)|,
       g_k = G prod_(j=k)^(K-1) (|be0 + j| + sqrt(2) R)
           >= sup |1 / Gamma(be + k)|,

   and T the bound from the terms for 1F1(al + K, be + K, z), where
   the real part of the lower parameter is at least 2 - R. */

/* log2(1 / R), the largest log2 |eps|, the largest order J, the largest
   number of terms up to the pole */
#define U_NEAR_RADIUS_EXP 2
#define U_NEAR_MAX_EXP (-12)
#define U_NEAR_MAX_ORDER 60
#define U_NEAR_MAX_POLE 20000

static void
u_near_series_majorant(mag_t S, const acb_t al, slong be0, const mag_t R, const acb_t z)
{
    slong K = FLINT_MAX(0, 2 - be0), k;
    mag_t zm, G, g, t, u, Rd;
    mag_ptr cs;
    acb_t v, w;

    if (K > U_NEAR_MAX_POLE || !acb_is_finite(al) || !acb_is_finite(z))
    {
        mag_inf(S);
        return;
    }

    mag_init(zm); mag_init(G); mag_init(g); mag_init(t); mag_init(u); mag_init(Rd);
    acb_init(v); acb_init(w);
    cs = _mag_vec_init(K + 1);

    acb_get_mag(zm, z);

    /* sqrt(2) R, rounded up */
    mag_set_d(Rd, 1.4143);
    mag_mul(Rd, Rd, R);

    /* G and T */
    acb_set_si(w, be0);
    acb_add_si(w, w, K, 64);
    acb_add_error_mag(w, R);
    acb_rgamma(v, w, 30);
    acb_get_mag(G, v);
    acb_add_si(v, al, K, 64);
    _acb_hypgeom_inexact_terms_bound(t, v, 1, w, 1, z, 30000, NULL);

    mag_one(cs);
    for (k = 0; k < K; k++)
    {
        acb_add_si(v, al, k, 64);
        acb_get_mag(u, v);
        mag_mul(cs + k + 1, cs + k, u);
        mag_mul(cs + k + 1, cs + k + 1, zm);
        mag_div_ui(cs + k + 1, cs + k + 1, k + 1);
    }

    mag_mul(S, cs + K, G);
    mag_mul(S, S, t);

    mag_set(g, G);
    for (k = K - 1; k >= 0; k--)
    {
        mag_set_ui(u, FLINT_ABS(be0 + k));
        mag_add(u, u, Rd);
        mag_mul(g, g, u);
        mag_addmul(S, cs + k, g);
    }

    mag_clear(zm); mag_clear(G); mag_clear(g); mag_clear(t); mag_clear(u); mag_clear(Rd);
    acb_clear(v); acb_clear(w);
    _mag_vec_clear(cs, K + 1);
}

/* M as above for the integer n */
static void
u_near_majorant(mag_t M, const acb_t a, slong n, const acb_t z)
{
    mag_t R, S, t;
    acb_t D, v, w;

    mag_init(R); mag_init(S); mag_init(t);
    acb_init(D); acb_init(v); acb_init(w);

    mag_set_ui_2exp_si(R, 1, -U_NEAR_RADIUS_EXP);
    acb_set_si(D, n);
    acb_add_error_mag(D, R);

    /* |A| <= |1 / Gamma(a - b + 1)| sup |1F~1(a, b, z)| */
    acb_sub(w, a, D, 64);
    acb_add_ui(w, w, 1, 64);
    acb_rgamma(v, w, 30);
    acb_get_mag(M, v);
    u_near_series_majorant(S, a, n, R, z);
    mag_mul(M, M, S);

    /* |B| <= |z^(1-b)| / |Gamma(a)| sup |1F~1(a - b + 1, 2 - b, z)| */
    acb_rgamma(v, a, 30);
    acb_get_mag(t, v);
    u_near_series_majorant(S, w, 2 - n, R, z);
    mag_mul(t, t, S);
    acb_sub_ui(v, D, 1, 64);
    acb_neg(v, v);
    acb_pow(v, z, v, 30);
    acb_get_mag(S, v);
    mag_mul(t, t, S);

    mag_add(M, M, t);

    /* pi sqrt(2) */
    mag_set_d(t, 4.4429);
    mag_mul(M, M, t);

    if (!acb_is_finite(v))
        mag_inf(M);

    mag_clear(R); mag_clear(S); mag_clear(t);
    acb_clear(D); acb_clear(v); acb_clear(w);
}

/* The order J of the expansion for b at the precision prec, or 0 if it
   is not used: b is an exact integer, or not within 2^-12 of one; a
   point or a ball without the integer for which the formula with more
   precision costs less; or more than 60 terms.  (On a ball of radius
   rad at the distance |eps| the formula gives log2(|eps| / rad) bits
   at best, whatever the precision is.) */
static slong
u_near_order(const acb_t b, slong prec)
{
    arf_srcptr m = arb_midref(acb_realref(b));
    acb_t eps;
    mag_t e;
    slong n, J = 0;
    double lq;

    if (!acb_is_finite(b) || acb_is_int(b) || arf_cmpabs_2exp_si(m, 40) >= 0)
        return 0;

    acb_init(eps);
    mag_init(e);

    n = arf_get_si(m, ARF_RND_NEAR);
    acb_sub_si(eps, b, n, FLINT_MAX(prec, acb_bits(b)) + 64);
    acb_get_mag(e, eps);

    if (mag_cmp_2exp_si(e, U_NEAR_MAX_EXP) < 0)
    {
        /* log2(1 / q), q = |eps| / R */
        lq = -mag_get_d_log2_approx(e) - U_NEAR_RADIUS_EXP - 0.01;
        J = (slong) ceil((prec + 16.0) / lq) + 1;

        /* (A point, or a ball which does not contain the integer: the
           formula at the midpoint with the bits which it loses, and a
           bound for the variation, cost less than a long expansion.
           The orders are where the timings of the two cross at 64,
           333 and 3333 bits.  A ball at low precision on which the
           formula loses accuracy is the exception: there the function
           accepts the result of the formula if it has lost less than
           30 bits, which the expansion up to the order 12 does not
           lose.) */
        if (J > U_NEAR_MAX_ORDER)
        {
            J = 0;
        }
        else
        {
            /* (a ball without the integer on which the formula keeps
               the accuracy, log2(|eps| / rad) bits, is like a point) */
            int point = acb_is_exact(b), contains = 0;

            if (!point)
            {
                mag_t lo, r;

                mag_init(lo);
                mag_init(r);
                acb_get_mag_lower(lo, eps);
                mag_hypot(r, arb_radref(acb_realref(b)), arb_radref(acb_imagref(b)));
                contains = mag_is_zero(lo);
                if (!contains)
                {
                    mag_mul_2exp_si(r, r, prec - 4);
                    point = (mag_cmp(r, lo) <= 0);
                }
                mag_clear(lo);
                mag_clear(r);
            }

            if (point)
            {
                if (J > ((prec < 200) ? 2 : (prec < 1500) ? 3 : 6))
                    J = 0;
            }
            else if (!contains && J > ((prec < 200) ? 12 : (prec < 1500) ? 4 : 6))
            {
                J = 0;
            }
        }
    }

    acb_clear(eps);
    mag_clear(e);

    return J;
}

/* U(a, b, z) by the expansion in b at the nearest integer; returns 0 if
   it does not apply */
static int
u_near_integer(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_poly_t aa, bb, zz, s;
    acb_t eps, r;
    mag_t e, q, M, rem;
    slong n, J, wp;
    double lq;
    int real, iter, ok = 0;

    J = u_near_order(b, prec);

    if (J == 0 || !acb_is_finite(a) || !acb_is_finite(z) || acb_contains_zero(z))
        return 0;

    acb_poly_init(aa); acb_poly_init(bb); acb_poly_init(zz); acb_poly_init(s);
    acb_init(eps); acb_init(r);
    mag_init(e); mag_init(q); mag_init(M); mag_init(rem);

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z) && arb_is_positive(acb_realref(z));

    n = arf_get_si(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    acb_sub_si(eps, b, n, FLINT_MAX(prec, acb_bits(b)) + 64);
    acb_get_mag(e, eps);
    mag_mul_2exp_si(q, e, U_NEAR_RADIUS_EXP);
    lq = -mag_get_d_log2_approx(q) - 0.01;

    u_near_majorant(M, a, n, z);

    if (mag_is_finite(M))
    {
        acb_poly_set_acb(aa, a);
        acb_poly_set_coeff_si(bb, 0, n);
        acb_poly_set_coeff_si(bb, 1, 1);
        acb_poly_set_acb(zz, z);

        for (iter = 0; iter < 3; iter++)
        {
            slong acc;
            double miss;

            wp = prec + 10 + 2 * J;
            acb_hypgeom_u_1f1_series(s, aa, bb, zz, J, wp);
            acb_poly_evaluate(r, s, eps, wp);

            /* M q^J / (1 - q) */
            mag_geom_series(rem, q, J);
            mag_mul(rem, rem, M);

            if (real)
            {
                arb_zero(acb_imagref(r));
                arb_add_error_mag(acb_realref(r), rem);
            }
            else
            {
                acb_add_error_mag(r, rem);
            }

            if (!acb_is_finite(r))
                break;

            ok = 1;
            acb_set(res, r);
            acc = acb_rel_accuracy_bits(r);

            if (acc >= prec - 8)
                break;

            /* (M is larger than the function by the cancellation of A
               and B: more terms, if the remainder is what is missing) */
            acb_get_mag_lower(e, r);
            if (mag_is_zero(e))
                break;
            miss = mag_get_d_log2_approx(rem) - mag_get_d_log2_approx(e) + prec;
            if (miss < 2.0)
                break;
            J += (slong) ceil((miss + 8.0) / lq);
            if (J > U_NEAR_MAX_ORDER)
                break;
        }
    }

    acb_poly_clear(aa); acb_poly_clear(bb); acb_poly_clear(zz); acb_poly_clear(s);
    acb_clear(eps); acb_clear(r);
    mag_clear(e); mag_clear(q); mag_clear(M); mag_clear(rem);

    return ok;
}

/* (n_max: the largest number of terms of the two series for
   noninteger b; by default, a multiple of the precision) */
static void
_acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec, slong n_max)
{
    n_max = FLINT_MAX(n_max, 50 + 10 * prec);

    /* b close to an integer: the expansion in b at the integer */
    if (u_near_integer(res, a, b, z, prec))
        return;

    if (acb_is_int(b))
    {
        acb_poly_t aa, bb, zz;

        acb_poly_init(aa);
        acb_poly_init(bb);
        acb_poly_init(zz);

        acb_poly_set_acb(aa, a);
        acb_poly_set_coeff_acb(bb, 0, b);
        acb_poly_set_coeff_si(bb, 1, 1);
        acb_poly_set_acb(zz, z);

        acb_hypgeom_u_1f1_series(zz, aa, bb, zz, 1, prec);
        acb_poly_get_coeff_acb(res, zz, 0);

        acb_poly_clear(aa);
        acb_poly_clear(bb);
        acb_poly_clear(zz);
    }
    else
    {
        acb_t t, u, v;
        acb_struct aa[3];

        acb_init(t);
        acb_init(u);
        acb_init(v);
        acb_init(aa + 0);
        acb_init(aa + 1);
        acb_init(aa + 2);

        acb_set(aa, a);
        acb_set(aa + 1, b);
        acb_one(aa + 2);
        _acb_hypgeom_pfq_direct_max(u, aa, 1, aa + 1, 2, z, prec, n_max);

        acb_sub(aa, a, b, prec);
        acb_add_ui(aa, aa, 1, prec);
        acb_sub_ui(aa + 1, b, 2, prec);
        acb_neg(aa + 1, aa + 1);
        _acb_hypgeom_pfq_direct_max(v, aa, 1, aa + 1, 2, z, prec, n_max);

        acb_sub_ui(aa + 1, b, 1, prec);

        /* rgamma(a-b+1) * gamma(1-b) * u */
        acb_rgamma(t, aa, prec);
        acb_mul(u, u, t, prec);
        acb_neg(t, aa + 1);
        acb_gamma(t, t, prec);
        acb_mul(u, u, t, prec);

        /* rgamma(a) * gamma(b-1) * z^(1-b) * v */
        acb_rgamma(t, a, prec);
        acb_mul(v, v, t, prec);
        acb_gamma(t, aa + 1, prec);
        acb_mul(v, v, t, prec);
        acb_neg(t, aa + 1);
        acb_pow(t, z, t, prec);
        acb_mul(v, v, t, prec);

        acb_add(res, u, v, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        acb_clear(aa + 0);
        acb_clear(aa + 1);
        acb_clear(aa + 2);
    }
}

void
acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    _acb_hypgeom_u_1f1(res, a, b, z, prec, 0);
}

static void
acb_hypgeom_u_choose(int * asymp, slong * wp,
    const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    double x, y, t, cancellation;
    double input_accuracy, direct_accuracy, asymp_accuracy;

    *asymp = 0;
    *wp = prec;

    input_accuracy = acb_rel_one_accuracy_bits(z);
    t = acb_rel_one_accuracy_bits(a);
    input_accuracy = FLINT_MIN(input_accuracy, t);
    t = acb_rel_one_accuracy_bits(b);
    input_accuracy = FLINT_MIN(input_accuracy, t);
    input_accuracy = FLINT_MAX(input_accuracy, 0.0);

    /* From here we ignore the values of a, b. Taking them into account is
       a possible future improvement... */

    /* Tiny |z|. */
    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 2) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 2) < 0))
    {
        *asymp = 0;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    /* Huge |z|. */
    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    asymp_accuracy = sqrt(x * x + y * y) * 1.44269504088896;

    /* The Kummer transformation gives less cancellation with the 1F1 series.
    if (x < 0.0)
    {
        *kummer = 1;
        x = -x;
    } */

    if (asymp_accuracy >= prec)
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    /* Assume U ~ 1, M ~ exp(|z|)  (there is cancellation both in the
       evaluation of M and in the linear combination) -- a better estimate
       would account for a, b. */
    cancellation = sqrt(x * x + y * y) * 1.44269504088896 + 5;

    direct_accuracy = input_accuracy - cancellation;

    if (direct_accuracy > asymp_accuracy)
    {
        *asymp = 0;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec + cancellation));
    }
    else
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
    }
}

void
_acb_hypgeom_u_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t t;
    arf_srcptr av, tv;

    av = arb_midref(acb_realref(a));

    /* Handle small polynomial cases without divisions. */
    /* todo: should incorporate a -> 1+a-b transformation, also... */
    if (acb_is_int(a) && arf_sgn(av) <= 0)
    {
        if (arf_cmpabs_ui(av, 30) < 0 ||
            (arf_cmpabs_ui(av, prec) < 0 &&
                ((acb_bits(b) < 0.1 * prec && acb_bits(z) < 0.1 * prec)
                    || acb_contains_zero(z))))
        {
            acb_hypgeom_u_si_rec(res, arf_get_si(av, ARF_RND_DOWN), b, z, prec);
            return;
        }
    }

    acb_init(t);
    acb_sub(t, a, b, prec);
    acb_add_ui(t, t, 1, prec);
    tv = arb_midref(acb_realref(t));

    /* todo: combine these conditions with the code below */
    if ((acb_is_int(a) && arf_sgn(av) <= 0) ||
        (acb_is_int(t) && arf_sgn(tv) <= 0) ||
        acb_hypgeom_u_use_asymp(z, prec))
    {
        acb_neg(t, a);
        acb_pow(t, z, t, prec);
        acb_hypgeom_u_asymp(res, a, b, z, -1, prec);
        acb_mul(res, res, t, prec);
    }
    else
    {
        slong wp;
        int asymp;

        acb_hypgeom_u_choose(&asymp, &wp, a, b, z, prec);

        if (asymp)
        {
            acb_neg(t, a);
            acb_pow(t, z, t, prec);
            acb_hypgeom_u_asymp(res, a, b, z, -1, wp);
            acb_mul(res, res, t, prec);
        }
        else
        {
            acb_hypgeom_u_1f1(res, a, b, z, wp);
            acb_set_round(res, res, prec);
        }
    }

    acb_clear(t);
}

/* ------------------------------------------------------------------------- */
/* Choice of algorithm from estimates of the costs                            */
/* ------------------------------------------------------------------------- */

/* See the comments for acb_hypgeom_m; the estimates are in estimates.c. */

#define U_INVLOG2 1.4426950408889634

typedef _acb_hypgeom_cd u_cd;

/* Numerical integration, replacing best if more accurate (the complex
   integral if the real one does not apply, or is less accurate than
   enough). */
static void
u_integrate(acb_t best, slong * bestacc, const acb_t a, const acb_t b, const acb_t z,
    slong iprec, slong enough)
{
    acb_t t;
    slong acc, racc = -WORD_MAX;

    acb_init(t);

    /* real arguments: on the real line (and for large parameters along
       paths of steepest descent) */
    if (acb_is_real(a) && acb_is_real(b) && acb_is_real(z) && arb_is_positive(acb_realref(z)))
    {
        arb_hypgeom_u_integration(acb_realref(t),
                acb_realref(a), acb_realref(b), acb_realref(z), iprec);
        arb_zero(acb_imagref(t));
        acc = racc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    /* complex arguments (or real ones outside of the domain of the real
       integral): along paths of steepest descent */
    if (racc < enough)
    {
        acb_hypgeom_u_integration(t, a, b, z, iprec);
        acc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    acb_clear(t);
}

/* U by the 1F1 series (asymp = 0) or the asymptotic series at precision
   wp */
static void
u_series(acb_t res, const acb_t a, const acb_t b, const acb_t z, int asymp, slong wp, slong n_max)
{
    if (asymp)
    {
        acb_t t;
        acb_init(t);
        acb_neg(t, a);
        acb_pow(t, z, t, wp);
        acb_hypgeom_u_asymp(res, a, b, z, -1, wp);
        acb_mul(res, res, t, wp);
        acb_clear(t);
    }
    else
    {
        _acb_hypgeom_u_1f1(res, a, b, z, wp, n_max);
    }
}

/* ------------------------------------------------------------------------- */
/* Inexact input                                                              */
/* ------------------------------------------------------------------------- */

/* Bound for |U| on balls A, B, Z with Re(A) > 0 and Re(Z) > 0, from

       Gamma(a) U(a, b, z) = int_0^inf e^(-z t) t^(a-1) (1 + t)^(b-a-1) dt:

   with the real parts al, be, x of a, b, z,

       |U(a, b, z)| <= Gamma(al) / |Gamma(a)| U(al, be, x),

   Gamma(al) / |Gamma(al + i y)| <= exp(y^2 (1 / al + 1 / al^2) / 2) (see
   the bound for 1F1), and the integral Gamma(al) U(al, be, x) decreases
   with al (the integrand is (t / (1 + t))^al times a function of be)
   and with x, and increases with be: its largest value on the box of
   real parts is at a corner.  This bound is close to |U| unless z or
   the parameters are far from real. */
static void
u_majorant_euler(mag_t M, const acb_t A, const acb_t B, const acb_t Z)
{
    acb_t t, al, be, x;
    arb_t g, u, w;
    arf_t c;
    slong prec = 40;

    acb_init(t); acb_init(al); acb_init(be); acb_init(x);
    arb_init(g); arb_init(u); arb_init(w);
    arf_init(c);

    if (arb_is_positive(acb_realref(A)) && arb_is_positive(acb_realref(Z))
        && acb_is_finite(A) && acb_is_finite(B) && acb_is_finite(Z))
    {
        arb_get_lbound_arf(c, acb_realref(A), prec);
        acb_zero(al); arb_set_arf(acb_realref(al), c);
        arb_get_ubound_arf(c, acb_realref(B), prec);
        acb_zero(be); arb_set_arf(acb_realref(be), c);
        arb_get_lbound_arf(c, acb_realref(Z), prec);
        acb_zero(x); arb_set_arf(acb_realref(x), c);

        /* exp(y^2 (1 / al + 1 / al^2) / 2) */
        arb_inv(u, acb_realref(al), prec);
        arb_mul(w, u, u, prec);
        arb_add(u, u, w, prec);
        arb_mul(w, acb_imagref(A), acb_imagref(A), prec);
        arb_mul(u, u, w, prec);
        arb_mul_2exp_si(u, u, -1);
        arb_exp(g, u, prec);

        /* 1 / Gamma(al) on the interval, times the integral at the corner */
        arb_rgamma(u, acb_realref(A), prec);
        arb_mul(g, g, u, prec);
        arb_gamma(u, acb_realref(al), prec);
        arb_mul(g, g, u, prec);

        acb_hypgeom_u(t, al, be, x, prec);
        acb_mul_arb(t, t, g, prec);
        acb_get_mag(M, t);
    }
    else
    {
        mag_inf(M);
    }

    acb_clear(t); acb_clear(al); acb_clear(be); acb_clear(x);
    arb_clear(g); arb_clear(u); arb_clear(w);
    arf_clear(c);
}

static void
u_eval_exact(acb_t res, acb_srcptr w, slong prec, void * FLINT_UNUSED(ctx))
{
    acb_hypgeom_u(res, w, w + 1, w + 2, prec);
}

/* The bounds for |U| on balls (see _acb_hypgeom_inexact) */
static void
u_majorant(mag_t M, acb_srcptr w, int stage, void * FLINT_UNUSED(ctx))
{
    acb_t t;
    acb_init(t);

    if (stage == 0)
    {
        /* from the integral, also for z^(1-b) U(a - b + 1, 2 - b, z) */
        acb_t a2, b2;
        mag_t M2;

        acb_init(a2);
        acb_init(b2);
        mag_init(M2);

        acb_sub(a2, w, w + 1, 40);
        acb_add_ui(a2, a2, 1, 40);
        acb_sub_ui(b2, w + 1, 2, 40);
        acb_neg(b2, b2);
        u_majorant_euler(M2, a2, b2, w + 2);

        if (mag_is_finite(M2))
        {
            acb_sub_ui(b2, w + 1, 1, 40);
            acb_neg(b2, b2);
            acb_pow(t, w + 2, b2, 40);
            acb_get_mag(M, t);
            mag_mul(M2, M2, M);
        }

        u_majorant_euler(M, w, w + 1, w + 2);
        mag_min(M, M, M2);

        acb_clear(a2);
        acb_clear(b2);
        mag_clear(M2);
    }
    else if (stage == 1)
    {
        /* the series on the balls */
        _acb_hypgeom_u_1f1(t, w, w + 1, w + 2, 64, 0);
        acb_get_mag(M, t);

        /* (a ball for b which contains an integer n, inside of the
           disc of radius 1/4 around it: the bound from the terms of
           the limit) */
        if (!mag_is_finite(M) && acb_is_finite(w + 1)
            && arf_cmpabs_2exp_si(arb_midref(acb_realref(w + 1)), 40) < 0)
        {
            slong n = arf_get_si(arb_midref(acb_realref(w + 1)), ARF_RND_NEAR);
            mag_t r;

            mag_init(r);
            acb_sub_si(t, w + 1, n, 64);
            acb_get_mag(r, t);
            if (mag_cmp_2exp_si(r, -U_NEAR_RADIUS_EXP) <= 0)
                u_near_majorant(M, w, n, w + 2);
            mag_clear(r);
        }
    }
    else if (stage == 2)
    {
        /* the asymptotic series on the balls */
        mag_t r;
        mag_init(r);
        acb_get_mag_lower(r, w + 2);
        if (mag_cmp_2exp_si(r, 3) >= 0)
        {
            acb_t u;
            acb_init(u);
            acb_neg(u, w);
            acb_pow(u, w + 2, u, 64);
            acb_hypgeom_u_asymp(t, w, w + 1, w + 2, -1, 64);
            acb_mul(t, t, u, 64);
            acb_get_mag(M, t);
            acb_clear(u);
        }
        else
            mag_inf(M);
        mag_clear(r);
    }
    else
    {
        /* the integrals along paths of steepest descent */
        acb_hypgeom_u_integration(t, w, w + 1, w + 2, 30);
        acb_get_mag(M, t);
    }

    acb_clear(t);
}

/* The function on inexact input from its value at the midpoints and a
   bound for the variation (inexact.c).  Returns 0 if there is no
   result, 2 if the enclosure is close to the range of the function. */
static int
u_inexact(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec, slong goal,
    double cancel, double exp_ratio, int lazy)
{
    acb_struct args[3];
    double R[3], sa, sb, sz, spar, na, nb, nz, az, dist, q;
    int n, real;

    args[0] = *a;
    args[1] = *b;
    args[2] = *z;

    n = !acb_is_exact(a) + !acb_is_exact(b) + !acb_is_exact(z);

    na = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    sa = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    nb = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    sb = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
    nz = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    sz = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    az = hypot(nz, sz);

    if (!(hypot(na, sa) < 1e15) || !(hypot(nb, sb) < 1e15) || !(az < 1e15) || !(az > 1e-15))
        return 0;

    /* The radii of the polydisc: about the scales on which the function
       varies by a factor e.  The logarithmic derivatives in a and b are
       of the size of log |parameters| (or log |z|), and in z of size 1
       or |a / z|, |b / z|.  The disc for z must not meet the cut. */
    spar = 2.0 + log(2.0 + hypot(na, sa) + hypot(nb, sb) + az + 1.0 / az);
    R[0] = R[1] = 1.0 / (n * spar);
    R[2] = 1.0 / (n * (1.0 + (hypot(na, sa) + hypot(nb, sb)) / az));
    dist = (nz >= 0.0) ? az : fabs(sz);
    R[2] = FLINT_MIN(R[2], 0.05 * dist);

    /* (inside of the domains of the integrals, if the midpoint is) */
    if (na > 0.0)
        R[0] = FLINT_MIN(R[0], 0.05 * na);
    q = na - nb + 1.0;
    if (q > 0.0)
    {
        R[0] = FLINT_MIN(R[0], 0.025 * q);
        R[1] = FLINT_MIN(R[1], 0.025 * q);
    }
    q = nb - na;
    if (q > 0.0)
    {
        R[0] = FLINT_MIN(R[0], 0.025 * q);
        R[1] = FLINT_MIN(R[1], 0.025 * q);
    }

    /* (The terms of the series have poles at integers b, which cancel.
       For the bound from the series on balls, the disc for b does not
       contain an integer; for the bound from the terms of the limit,
       its square is inside of the disc of radius 1/4 around the
       integer.) */
    q = hypot(nb - floor(nb + 0.5), sb);
    if (fabs(sb) < R[1])
        R[1] = FLINT_MIN(R[1], FLINT_MAX(0.1 * fabs(nb - floor(nb + 0.5)), (0.249 - q) / 1.4143));

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z) && arb_is_positive(acb_realref(z));

    return _acb_hypgeom_inexact(res, args, 3, R, u_eval_exact, u_majorant, 2,
                NULL, real, prec, goal, cancel, exp_ratio, lazy, 0);
}

#define U_SERIES 0
#define U_ASYMP 1
#define U_INTEGRATION 2

static int
u_is_nonpositive_int(const acb_t x)
{
    return acb_is_int(x) && arf_sgn(arb_midref(acb_realref(x))) <= 0;
}

/* The Taylor coefficients in z for _acb_hypgeom_inexact_taylor:
   (-1)^j (a)_j / j! U(a + j, b + j, z). */
typedef struct
{
    const acb_struct * a;
    const acb_struct * b;
}
u_taylor_ctx;

#define U_TAYLOR_MAX_ORDER 60
#define U_TAYLOR_MAX_VALUES 8

static void
u_taylor_coeff(acb_t res, slong j, const acb_t z, slong prec, int crude, void * vctx)
{
    u_taylor_ctx * ctx = (u_taylor_ctx *) vctx;
    acb_t a, b, t;
    arb_t f;
    slong wp = prec + 10;

    acb_init(a);
    acb_init(b);
    acb_init(t);
    arb_init(f);

    acb_rising_ui(t, ctx->a, j, wp);

    if (acb_is_zero(t))
    {
        /* (a polynomial of degree less than j) */
        acb_zero(res);
    }
    else
    {
        arb_fac_ui(f, j, wp);
        acb_div_arb(t, t, f, wp);
        if (j % 2)
            acb_neg(t, t);

        /* (exact, unless the parameters have more bits than that) */
        acb_add_ui(a, ctx->a, j, wp + 200);
        acb_add_ui(b, ctx->b, j, wp + 200);

        if (!acb_is_exact(a) || !acb_is_exact(b))
            acb_indeterminate(res);
        else if (crude)
            _acb_hypgeom_u_nointegration(res, a, b, z, prec);
        else
            acb_hypgeom_u(res, a, b, z, prec);

        acb_mul(res, res, t, wp);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
    arb_clear(f);
}


/* c_{j+2} from c_j and c_{j+1} at m, for a solution of Kummer's
   equation z w'' + (b - z) w' - a w = 0:
   m (j+2)(j+1) c_{j+2} = (j + a) c_j - (j+1)(j + b - m) c_{j+1}. */
static void
u_taylor_recur(acb_t res, const acb_t c0, const acb_t c1, slong j, const acb_t m,
    slong prec, void * vctx)
{
    u_taylor_ctx * ctx = (u_taylor_ctx *) vctx;
    acb_t s, t;

    acb_init(s);
    acb_init(t);

    acb_add_si(s, ctx->a, j, prec);
    acb_mul(s, s, c0, prec);

    acb_add_si(t, ctx->b, j, prec);
    acb_sub(t, t, m, prec);
    acb_mul(t, t, c1, prec);
    acb_mul_si(t, t, j + 1, prec);
    acb_sub(s, s, t, prec);

    acb_mul_si(t, m, (j + 2) * (j + 1), prec);
    acb_div(res, s, t, prec);

    acb_clear(s);
    acb_clear(t);
}

/* Whether U with a or a - b + 1 a nonpositive integer is analytic on
   the box z: a polynomial in z for such a, and z^(1 - b) times one
   otherwise, which is rational for integer b and has the cut
   (-inf, 0] otherwise. */
static int
u_poly_analytic(const acb_t a, const acb_t b, const acb_t z)
{
    if (u_is_nonpositive_int(a))
        return 1;

    if (acb_contains_zero(z))
        return 0;

    if (arb_is_positive(acb_realref(z)) || !arb_contains_zero(acb_imagref(z)))
        return 1;

    return acb_is_int(b);
}

/* The number of terms of the polynomial in 1/z (a or a - b + 1 a
   nonpositive integer; the smaller number of terms), or 0. */
static double
u_poly_terms(const acb_t a, const acb_t b, slong prec)
{
    double n = 0.0, m;
    acb_t t;

    acb_init(t);

    if (u_is_nonpositive_int(a))
        n = 1.0 - arf_get_d(arb_midref(acb_realref(a)), ARF_RND_UP);

    acb_sub(t, a, b, prec);
    acb_add_ui(t, t, 1, prec);
    if (acb_is_exact(b) && u_is_nonpositive_int(t))
    {
        m = 1.0 - arf_get_d(arb_midref(acb_realref(t)), ARF_RND_UP);
        if (n == 0.0 || m < n)
            n = m;
    }

    acb_clear(t);
    return n;
}

/* (the cost of nterms terms at the precision wp, in microseconds) */
static double
u_poly_cost(const acb_t a, const acb_t b, const acb_t z, double wp, double nterms, int real)
{
    double pbits, zterm;

    pbits = acb_bits(a);
    pbits = FLINT_MAX(pbits, acb_bits(b));

    /* (an argument of many bits: see f_poly_cost in 2f1.c) */
    zterm = (acb_bits(z) > 64) ? (real ? 1.3e-4 : 4e-4) * wp : 0.0;

    return nterms * (_acb_hypgeom_est_series_term_cost(wp, pbits, real) + zterm)
        + 40.0 + 0.01 * pow(wp, 1.5);
}

/* A polynomial with exact parameters and argument whose terms cancel
   by more bits than the precision provides: the same computation with
   the missing bits, or with twice the precision when the result is all
   noise, while the cost stays within the budget (in microseconds: 50 ms,
   or half the estimated cost of the integral which follows if the
   polynomial fails, and if it has been tried, 16 times the cost of the
   polynomial at the precision of the caller, as for 2F1).  The precision is limited as well,
   since a polynomial which vanishes at z is not accurate at any
   precision. */
static void
u_poly_refine(acb_t best, slong * bestacc, const acb_t a, const acb_t b, const acb_t z,
    slong prec, slong accept, double nterms, int real, double maxwp, double budget)
{
    acb_t t;
    slong wp = prec, wp2, acc = *bestacc;
    double limit, spent = 0.0, cost;
    int iter;

    acb_init(t);

    limit = 16.0 * prec + 8.0 * nterms + 1000.0;
    limit = FLINT_MIN(limit, maxwp);

    for (iter = 0; iter < 30 && *bestacc < accept; iter++)
    {
        wp2 = (acc >= 4) ? wp + (accept - acc) + 10 + wp / 64 : 2 * wp;

        cost = u_poly_cost(a, b, z, wp2, nterms, real);

        if (wp2 > limit || spent + cost > budget)
            break;

        spent += cost;
        wp = wp2;
        _acb_hypgeom_u_nointegration(t, a, b, z, wp);
        acc = acb_rel_accuracy_bits(t);

        if (acc > *bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            *bestacc = acc;
        }
    }

    acb_clear(t);
}

void
acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t best, t;
    u_cd A, B, Z, A2, B2, lz;
    double L[2], P[2], pref[2], kpk[2], Lmax, lu, lu2, pbits, c_asymp, n_asymp, asymp_wp;
    double budget, spent, zabs, cref = 1e300, spare;
    slong goal, accept, acc, bestacc, in_acc, in_min;
    int real, exact, int_ok, tried_asymp, tried_int, series_runs, iter, have_model, loop_first = 0;
    int lazy_input;
    int polynomial, intb, moderate;
    slong near;
    double limit_factor;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(z))
    {
        _acb_hypgeom_u_nointegration(res, a, b, z, prec);
        return;
    }

    acb_init(best);
    acb_init(t);

    real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z) && arb_is_positive(acb_realref(z));
    exact = acb_is_exact(a) && acb_is_exact(b) && acb_is_exact(z);
    lazy_input = !exact && !(_acb_hypgeom_inexact_is_coarse(a, prec)
        || _acb_hypgeom_inexact_is_coarse(b, prec) || _acb_hypgeom_inexact_is_coarse(z, prec));

    /* The accuracy that the input is good for: the relative condition
       number for an inexact z, a or b is about its absolute value. */
    in_min = _acb_hypgeom_est_input_accuracy(z);
    in_acc = _acb_hypgeom_est_input_accuracy(a);
    in_min = FLINT_MIN(in_min, in_acc);
    in_acc = _acb_hypgeom_est_input_accuracy(b);
    in_min = FLINT_MIN(in_min, in_acc);
    goal = FLINT_MIN(prec, in_min);
    goal = FLINT_MAX(goal, 2);

    /* The bits by which inexact input is more accurate than the
       precision: series which cancel by fewer bits are evaluated on
       the input with more precision, as for exact input. */
    spare = exact ? 1e300 : (double) in_min - (double) prec;

    /* (such input is not the result of rounding at this precision:
       there is no reason to expect a caller to provide more) */
    if (spare >= 30.0)
        lazy_input = 0;

    /* Nearly full accuracy, or what the input is good for. */
    accept = goal - (exact ? 10 : 12);
    bestacc = -WORD_MAX;
    acb_indeterminate(best);

    A.re = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    A.im = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    B.re = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
    B.im = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
    Z.re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    Z.im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    zabs = hypot(Z.re, Z.im);

    /* polynomials in 1/z (times a power of z) */
    acb_sub(t, a, b, prec);
    acb_add_ui(t, t, 1, prec);
    polynomial = u_is_nonpositive_int(a) || u_is_nonpositive_int(t);

    /* (an integer b, or b so close to an integer that the series are
       expanded in b: the limit with power series, whose cost per term
       grows with the order of the expansion) */
    near = u_near_order(b, FLINT_MIN(prec, goal + 20));
    intb = acb_is_int(b) || near > 0;
    limit_factor = (real ? 3.5 : 1.4) * ((near > 0) ? 0.6 + 0.35 * near : 1.0);

    have_model = hypot(A.re, A.im) < 1e12 && hypot(B.re, B.im) < 1e12 && zabs < 1e12
        && zabs > 1e-12 && !polynomial;

    /* A parameter a or a - b + 1 of the asymptotic series of moderate
       size: the choice by the size of z alone is usually right. */
    moderate = hypot(A.re, A.im) <= 30.0 || hypot(A.re - B.re + 1.0, A.im - B.im) <= 30.0;

    /* (Not if the other one is large, that choice is the asymptotic
       series and its terms do not get small enough, as for
       U(1, 1000 + 1/2, 1000): the attempt would cost several times as
       much as the series.) */
    if (have_model && moderate && (hypot(A.re, A.im) > 30.0
            || hypot(A.re - B.re + 1.0, A.im - B.im) > 30.0))
    {
        int asymp = acb_hypgeom_u_use_asymp(z, prec);

        if (!asymp)
        {
            slong cwp;
            acb_hypgeom_u_choose(&asymp, &cwp, a, b, z, prec);
        }

        if (asymp)
        {
            double n0, more, cancel;

            if (!_acb_hypgeom_est_u_asymp(&n0, &more, &cancel, A, B, Z, (goal + 3.0) / U_INVLOG2, real))
                moderate = 0;
        }
    }

    if (!have_model || moderate)
    {
        _acb_hypgeom_u_nointegration(best, a, b, z, prec);
        bestacc = acb_rel_accuracy_bits(best);

        /* (for inexact input that is as accurate as the precision, a
           few bits less are accepted up to 1000 bits: the function at
           the midpoints with a bound for the variation costs several
           times as much there) */
        if (bestacc >= accept - ((lazy_input && prec <= 1000) ? 18 : 0) || !have_model)
        {
            if (bestacc < accept && goal > 2)
            {
                double nterms = polynomial ? u_poly_terms(a, b, prec) : 0.0;
                slong iprec = (goal < prec - 30) ? goal + 30 : prec;

                if (nterms > 0.0)
                {
                    /* A polynomial which cancels: an integral, if it
                       costs less than a few attempts with more
                       precision (it returns at once where it does not
                       apply); with exact parameters and argument, the
                       polynomial with more precision, and otherwise
                       the polynomial at the midpoints with a bound for
                       the variation (see inexact.c); the integral if
                       that has failed. */
                    /* (the Euler integral for Re(a) >= 1, or two
                       integrals along complex paths) */
                    double cint = (A.re >= 1.0) ? 1.3 * _acb_hypgeom_est_integration_cost(iprec, real)
                        : (real ? 1.3 : 2.6) * _acb_hypgeom_est_integration_cost(iprec, 0);
                    double cpoly = u_poly_cost(a, b, z, 2.0 * prec, nterms, real);
                    int tried = 0;

                    /* (input which is more accurate than the
                       precision: first the polynomial with more
                       precision, as far as the input allows) */
                    if (!exact && spare >= 30.0)
                        u_poly_refine(best, &bestacc, a, b, z, prec, accept, nterms, real,
                            prec + spare - 10.0, 5e4);

                    if (bestacc < goal / 2 && cint < 4.0 * cpoly)
                    {
                        u_integrate(best, &bestacc, a, b, z, iprec,
                            exact ? goal / 2 : FLINT_MIN(accept, goal - 20));
                        tried = 1;
                    }

                    /* Only the argument is inexact: the Taylor
                       expansion at its midpoint, whose terms are
                       polynomials of the same kind with exact
                       arguments. */
                    if (!exact && bestacc < accept && acb_is_exact(a) && acb_is_exact(b)
                        && u_poly_analytic(a, b, z))
                    {
                        u_taylor_ctx ctx;
                        slong order;

                        ctx.a = a;
                        ctx.b = b;

                        if (_acb_hypgeom_inexact_taylor(t, &order, z, u_taylor_coeff,
                                u_taylor_recur, &ctx, real, prec, goal,
                                U_TAYLOR_MAX_ORDER, U_TAYLOR_MAX_VALUES))
                        {
                            acc = acb_rel_accuracy_bits(t);

                            if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
                            {
                                acb_swap(best, t);
                                bestacc = acc;
                            }

                            acb_set_round(res, best, prec);
                            acb_clear(best);
                            acb_clear(t);
                            return;
                        }
                    }

                    if (exact)
                    {
                        if (bestacc < accept)
                            u_poly_refine(best, &bestacc, a, b, z, prec, accept, nterms, real, 1e300,
                                FLINT_MAX(5e4, tried ? 16.0 * u_poly_cost(a, b, z, prec, nterms, real)
                                    : 0.5 * cint));
                    }
                    else if (bestacc < accept - 18)
                    {
                        if (u_inexact(t, a, b, z, prec, goal, 0.0,
                                5.0 * _acb_hypgeom_est_integration_cost(64, 0) / cpoly, lazy_input) != 0)
                        {
                            acc = acb_rel_accuracy_bits(t);

                            if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
                            {
                                acb_swap(best, t);
                                bestacc = acc;
                            }
                        }
                    }

                    if (bestacc < goal / 2 && !tried)
                        u_integrate(best, &bestacc, a, b, z, iprec, goal / 2);
                }
                else if (bestacc < goal / 2 && !polynomial)
                {
                    /* (without estimates, integration if less than
                       half of the bits are accurate) */
                    u_integrate(best, &bestacc, a, b, z, iprec, goal / 2);
                }
            }

            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }
    }

    /* The two terms of
       U = Gamma(1-b) / Gamma(a-b+1) 1F1(a, b, z)
           + Gamma(b-1) / Gamma(a) z^(1-b) 1F1(a-b+1, 2-b, z):
       the largest terms of the series, with the factors (without their
       poles at integers b, which cancel). */
    A2.re = A.re - B.re + 1.0; A2.im = A.im - B.im;
    B2.re = 2.0 - B.re; B2.im = -B.im;
    lz = _acb_hypgeom_cd_log(Z);

    P[0] = _acb_hypgeom_est_1f1_peak(kpk, A, B, Z);
    P[1] = _acb_hypgeom_est_1f1_peak(kpk + 1, A2, B2, Z);
    pref[0] = _acb_hypgeom_d_lgamma_re(1.0 - B.re, -B.im, 1) - _acb_hypgeom_est_lgamma_re(A2);
    pref[1] = _acb_hypgeom_d_lgamma_re(B.re - 1.0, B.im, 1) - _acb_hypgeom_est_lgamma_re(A)
        + (1.0 - B.re) * lz.re + B.im * lz.im;
    L[0] = pref[0] + P[0];
    L[1] = pref[1] + P[1];
    Lmax = FLINT_MAX(L[0], L[1]);

    _acb_hypgeom_est_u_logabs(&lu, &lu2, A, B, Z, Lmax + 2.0);

    pbits = FLINT_MAX(arf_bits(arb_midref(acb_realref(a))), arf_bits(arb_midref(acb_imagref(a))));
    pbits = FLINT_MAX(pbits, arf_bits(arb_midref(acb_realref(b))));
    pbits = FLINT_MAX(pbits, arf_bits(arb_midref(acb_imagref(b))));
    pbits = FLINT_MAX(pbits, arf_bits(arb_midref(acb_realref(z))));
    pbits = FLINT_MAX(pbits, arf_bits(arb_midref(acb_imagref(z))));

    /* The asymptotic series */
    {
        double more, cancel, need, base;

        need = (goal + 3.0) / U_INVLOG2;
        base = FLINT_MIN(prec, goal + 20) + 15.0
            + 1.5 * log(1.0 + hypot(A.re, A.im) + hypot(B.re, B.im) + zabs) * U_INVLOG2;
        asymp_wp = base;
        c_asymp = 1e300;

        if (_acb_hypgeom_est_u_asymp(&n_asymp, &more, &cancel, A, B, Z, need, real))
        {
            asymp_wp = FLINT_MAX(asymp_wp, n_asymp / 9.5);
            asymp_wp = FLINT_MAX(asymp_wp, (need + more) * U_INVLOG2 / 10.5 + 4.0);
            if (cancel > 0.0)
                asymp_wp = FLINT_MAX(asymp_wp, base + cancel * U_INVLOG2);

            if (asymp_wp <= 1e8)
                c_asymp = 20.0 + n_asymp * 1.5 * _acb_hypgeom_est_series_term_cost(asymp_wp, pbits, real)
                    + 0.002 * pow(asymp_wp, 1.5);
        }
    }

    budget = 5e4 + 5.0 * _acb_hypgeom_est_integration_cost(FLINT_MIN(prec, goal + 30), real);
    spent = 0.0;
    tried_asymp = tried_int = series_runs = 0;

    /* Re(a) >= 1 or Re(a - b + 1) >= 1 */
    acb_sub_ui(t, a, 1, prec);
    int_ok = arb_is_nonnegative(acb_realref(t));
    if (!int_ok)
    {
        acb_sub(t, a, b, prec);
        int_ok = arb_is_nonnegative(acb_realref(t));
    }

    /* Otherwise, from two integrals for Re(b - a) >= 1 or Re(a) <= 0
       (for z^(1-b) U(a - b + 1, 2 - b, z)) and z in the right
       half-plane, if the function oscillates (int_ok = 2). */
    if (!int_ok && arb_is_positive(acb_realref(z)))
    {
        acb_sub(t, b, a, prec);
        acb_sub_ui(t, t, 1, prec);
        if (arb_is_nonnegative(acb_realref(t)) && _acb_hypgeom_est_2u_oscillatory(A, B, Z))
            int_ok = 2;
        else if (arb_is_nonpositive(acb_realref(a)) && _acb_hypgeom_est_2u_oscillatory(A2, B2, Z))
            int_ok = 2;
    }

    /* Inexact input: the error bounds of the series on balls are
       inflated by the cancellation, which more precision does not help.
       The function is evaluated at the midpoints, with a bound for its
       variation (see inexact.c for what is attempted at which cost).
       What follows after this is for the case that there is no result
       (z on the cut, a ball containing a zero of the function, or huge
       numbers). */
inexact_input:

    if (!exact)
    {
        double cmain, cint, cexp, c, tw, T0, T1;
        slong wp = FLINT_MIN(prec, goal + 30);
        int status, lazy, euler, nothing, done = 0;

        /* the cost of the value: the series, the asymptotic series or
           the integral */
        cmain = c_asymp;
        c = (Lmax - lu) * U_INVLOG2;
        c = FLINT_MIN(FLINT_MAX(c, 0.0), 1e9);
        tw = wp + c + 20.0 + (wp + c) / 64.0;
        T0 = _acb_hypgeom_est_1f1_terms(A, B, Z, kpk[0], Lmax - tw / U_INVLOG2 - pref[0]);
        T1 = _acb_hypgeom_est_1f1_terms(A2, B2, Z, kpk[1], Lmax - tw / U_INVLOG2 - pref[1]);
        if (intb && B.re <= 0.0)
            T0 = FLINT_MAX(T0, 2.0 - B.re);
        if (intb && B2.re <= 0.0)
            T1 = FLINT_MAX(T1, 2.0 - B2.re);
        if (tw < 1e8 && T0 < 1e9 && T1 < 1e9)
            cmain = FLINT_MIN(cmain, (T0 + T1) * _acb_hypgeom_est_series_term_cost(tw, pbits, real)
                        * (intb ? limit_factor : 1.0) + 40.0 + 0.01 * pow(tw, 1.5));

        cint = (int_ok == 0) ? 1e300
            : (int_ok == 1) ? _acb_hypgeom_est_integration_cost(wp, real)
            : (real ? 1.0 : 2.0) * _acb_hypgeom_est_integration_cost(wp, 0);

        /* (the reference for the cost of an integral below) */
        cref = cmain;

        /* At low precision, a series that loses little (or the
           asymptotic series): on the input, as for exact input, which
           saves the bound.  (If that is not accurate, this is reached
           once more.) */
        if (loop_first == 0 && wp <= 100 && (c <= 8.0 || c_asymp < 1e300) && cmain < cint)
        {
            loop_first = 1;
            goto methods;
        }

        /* Input which is more accurate than the precision by more
           bits than the series lose: on the input with more precision,
           as for exact input.  (If that is not accurate, this is
           reached once more.) */
        if (loop_first == 0 && spare >= 20.0 + ((c_asymp < 1e300) ? 0.0 : c) && cmain < cint)
        {
            loop_first = 1;
            goto methods;
        }

        /* the bound from the integrals on complex balls */
        cexp = 5.0 * _acb_hypgeom_est_integration_cost(64, 0);

        /* (less accurate than the precision: more precision does not
           improve the result) */
        lazy = lazy_input;

        /* The series cancel by more bits than the input has, and the
           bound is out of proportion: the series give no result that
           is worth its cost; the caller can increase the precision.
           (Not in the domain of the integral with a positive
           integrand, where there is a bound that costs little.) */
        euler = arb_is_positive(acb_realref(z))
            && (arb_is_positive(acb_realref(a)) || A2.re > 0.0);
        nothing = lazy && c > goal + 11.0 && c_asymp >= 1e300
            && cexp > 2.0 * FLINT_MIN(cmain, cint) && !euler;

        /* Where the value would be computed by the integral: the
           integral on the input, which encloses the function without
           the loss of the series.  Where the series on the input
           give nothing (no result at all, or one from the midpoints
           with bounds that are computed whatever they cost), also an
           integral of up to four times the cost of the series at the
           precision at which they succeed. */
        if (int_ok != 0 && cint <= ((nothing || (!lazy && c > goal + 11.0)) ? 4.0 : 1.0) * cmain)
        {
            u_integrate(best, &bestacc, a, b, z,
                (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));

            /* (the value and a bound would cost two more integrals) */
            done = (bestacc >= accept - 18);
        }

        if (nothing || done)
        {
            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }

        cmain = FLINT_MIN(cmain, cint);

        status = u_inexact(t, a, b, z, prec, goal, c, cexp / cmain, lazy);

        if (status != 0)
        {
            acc = acb_rel_accuracy_bits(t);

            if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
            {
                acb_swap(best, t);
                bestacc = acc;
            }

            /* a bound that is not close to the function: the integral
               on the input, if that is not out of proportion */
            if (status == 1 && bestacc < accept && int_ok && (!lazy || cint <= 1.5 * cmain))
                u_integrate(best, &bestacc, a, b, z,
                    (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));

            acb_set_round(res, best, prec);
            acb_clear(best);
            acb_clear(t);
            return;
        }
    }

methods:

    for (iter = 0; iter < 8 && bestacc < accept && loop_first != 2; iter++)
    {
        double cost[3], T[2], cbest, c, base, tw;
        slong wps, cancel, guard;
        int s, j, method;

        c = (Lmax - lu) * U_INVLOG2;
        c = FLINT_MIN(FLINT_MAX(c, 0.0), 1e9);
        cancel = (slong) c;

        /* no more precision than an inexact input is good for */
        base = FLINT_MIN(prec, goal + 20);
        guard = 20 + (slong) ((base + c) / 64.0);
        tw = base + c + guard;

        /* (for integer b, the limit is computed with power series, and
           the number of terms is limited by the precision) */
        for (j = 0; j < 6; j++)
        {
            for (s = 0; s < 2; s++)
            {
                T[s] = _acb_hypgeom_est_1f1_terms(s ? A2 : A, s ? B2 : B, Z, kpk[s],
                            Lmax - tw / U_INVLOG2 - pref[s]);

                /* (the series with the lower parameter -m is summed
                   past the pole at k = m, where the terms of the limit
                   continue) */
                if (intb && (s ? B2.re : B.re) <= 0.0)
                    T[s] = FLINT_MAX(T[s], 2.0 - (s ? B2.re : B.re));
            }

            if (!intb || FLINT_MAX(T[0], T[1]) <= 9.0 * tw)
                break;

            tw = FLINT_MAX(1.1 * tw, FLINT_MAX(T[0], T[1]) / 9.0);
        }

        if (j == 6)
            tw = 1e300;

        if (tw > 1e8 || T[0] > 1e9 || T[1] > 1e9 || series_runs >= 4
            || (!exact && spare < c + 20.0 && c > ((int_ok && !tried_int)
                    ? FLINT_MAX(12.0, goal / 8.0) : goal / 2.0 + 10.0)))
        {
            cost[U_SERIES] = 1e300;
            wps = 0;
        }
        else
        {
            wps = (slong) tw;
            cost[U_SERIES] = (T[0] + T[1]) * _acb_hypgeom_est_series_term_cost(tw, pbits, real)
                    * (intb ? limit_factor : 1.0)
                + 40.0 + 0.01 * pow(tw, 1.5);     /* (gamma functions, power) */
        }

        if (cost[U_SERIES] > budget - spent)
            cost[U_SERIES] = 1e300;

        cost[U_ASYMP] = (tried_asymp >= 3 || c_asymp > budget - spent) ? 1e300 : c_asymp;
        cost[U_INTEGRATION] = (tried_int || !int_ok) ? 1e300
                    : (int_ok == 1) ? _acb_hypgeom_est_integration_cost(FLINT_MIN(prec, goal + 30), real)
                    : (real ? 1.0 : 2.0) * _acb_hypgeom_est_integration_cost(FLINT_MIN(prec, goal + 30), 0);

        /* (for inexact input that is as accurate as the precision:
           no integral at many times the cost of the other methods,
           whose results are limited by the input) */
        if (iter == 0 && !(cref < 1e300))
            cref = FLINT_MIN(cost[U_SERIES], cost[U_ASYMP]);
        if (lazy_input && cref < 1e300 && cost[U_INTEGRATION] > 4.0 * cref)
            cost[U_INTEGRATION] = 1e300;

        method = 0;
        cbest = cost[0];
        for (s = 1; s < 3; s++)
        {
            if (cost[s] < cbest)
            {
                cbest = cost[s];
                method = s;
            }
        }

        if (cbest >= 1e300)
            break;

        spent += cbest;

        if (method == U_INTEGRATION)
        {
            u_integrate(best, &bestacc, a, b, z,
                (goal < prec - 30) ? goal + 30 : prec, FLINT_MIN(accept, goal - 20));
            tried_int = 1;

            /* (the integral encloses the function on inexact input
               without the loss of the series) */
            if (!exact && acb_is_finite(best))
                break;

            continue;
        }
        else if (method == U_ASYMP)
        {
            slong macc;

            u_series(t, a, b, z, 1, (slong) asymp_wp, 0);
            acc = macc = acb_rel_accuracy_bits(t);
            tried_asymp++;

            /* A near miss: once more with the missing bits.  No
               accurate bit (the factor of the error bound was
               underestimated): with twice the precision, twice. */
            if (macc < accept && macc > 4 && tried_asymp == 1 && exact)
            {
                asymp_wp += accept - macc + 10;
                c_asymp *= 1.0 + (accept - macc + 10) / asymp_wp;
            }
            else if (macc <= 4 && tried_asymp <= 2 && exact && asymp_wp < 5e7)
            {
                c_asymp *= _acb_hypgeom_est_series_term_cost(2.0 * asymp_wp, pbits, real)
                    / _acb_hypgeom_est_series_term_cost(asymp_wp, pbits, real);
                asymp_wp *= 2.0;
                budget = spent + (budget - spent) / 3.0;
            }
            else
                tried_asymp = 3;
        }
        else
        {
            slong macc;

            u_series(t, a, b, z, 0, wps, (slong) (1.25 * FLINT_MAX(T[0], T[1]) + 100.0));
            acc = macc = acb_rel_accuracy_bits(t);
            series_runs++;

            if (!exact && spare < cancel + 20.0)
            {
                /* (more precision does not improve the error bounds
                   for inexact parameters) */
                series_runs = 4;
            }
            else if (macc < accept)
            {
                /* The cancellation was underestimated.  If the result
                   is not all noise, this shows by how much. */
                double loss;

                if (macc >= 4)
                {
                    loss = (double) cancel + (accept - macc) + 10.0;
                    lu = FLINT_MIN(lu, Lmax - loss / U_INVLOG2);
                }
                else if (lu2 < lu - 20.0)
                {
                    /* the other saddle point */
                    lu = lu2;
                    lu2 = 1e300;
                    budget = spent + (budget - spent) / 3.0;
                }
                else
                {
                    /* twice the precision */
                    loss = cancel + wps;
                    lu = FLINT_MIN(lu, Lmax - loss / U_INVLOG2);
                    budget = spent + (budget - spent) / 3.0;
                }
            }
            else if (acc < accept)
            {
                series_runs = 4;
            }
        }

        if (acc > bestacc || (acb_is_finite(t) && !acb_is_finite(best)))
        {
            acb_swap(best, t);
            bestacc = acc;
        }
    }

    if (loop_first == 1 && bestacc < accept)
    {
        loop_first = 2;
        goto inexact_input;
    }

    if (bestacc == -WORD_MAX)
    {
        _acb_hypgeom_u_nointegration(best, a, b, z, prec);
        bestacc = acb_rel_accuracy_bits(best);
    }

    acb_set_round(res, best, prec);

    acb_clear(best);
    acb_clear(t);
}
