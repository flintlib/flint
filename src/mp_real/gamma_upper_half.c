/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_mat.h"
#include "mag.h"
#include "mp_real.h"
#include "impl.h"

/* Gamma(1/2, y) = sqrt(pi) erfc(sqrt(y)) for exact rational y = p/q > 0,
   given as limb arrays (internally an fmpq, for the leaves).

   ASYMPTOTIC SERIES.  Gamma(1/2, y) = e^-y y^(-1/2) (sum_(j<N) (-1)^j
   (1/2)_j y^-j + R_N), and for real y > 0 the remainder is bounded by
   the first omitted term in absolute value (DLMF 8.11(ii), with
   N >= a - 1 = -1/2).  The terms decrease while j + 1/2 < y, the
   smallest being about e^-y, so the series serves when y exceeds about
   W log 2 for W bits.  The partial sum is a hypergeometric sum
   (mp_real_hypgeom_sum, with Q = 2p and R(j) = q - 2qj for y = p/q).

   CONTINUED FRACTION.  Gamma(1/2, y) = e^-y y^(1/2) / x with Legendre's
   fraction x = y + (1/2)/(1 + 1/(y + (3/2)/(1 + 2/(y + ...)))), that
   is b0 + a1/(b1 + a2/(b2 + ...)) with b_k = y (k even), 1 (k odd) and
   a_k = k/2 (DLMF 8.9.2 with a = 1/2).  The equivalence transformation
   by 2q (k even) and 2 (k odd) makes the elements integers:
   b_k = 2p (even), 2 (odd), a_k = 2qk, and the value 2qx.  All the
   elements are positive, so consecutive convergents enclose the value
   and no tail bound is needed.  The convergents come from the product
   of the matrices [b_k 1; a_k 0] (transposed, as a recurrence for
   _mp_real_linrec).  N steps give about 2.4 sqrt(N y) bits (measured;
   the asymptotic rate e^(-4 sqrt(N y)) is not reached at the sizes
   used here), so the fraction is preferred for large y when the
   asymptotic series is not accurate enough.

   LOWER SERIES.  Gamma(1/2, y) = sqrt(pi) - gamma(1/2, y) with
   gamma(1/2, y) = 2 sqrt(y) e^-y sum_(j>=0) (2y)^j / (3 5 ... (2j+1)),
   a hypergeometric series with positive terms (mp_real_hypgeom_series),
   evaluated with about 1.44 y extra bits for the cancellation. */

static void
_y_ball(mp_real_t y, const fmpq_t a, slong n)
{
    mp_real_t t;
    mp_real_init(t);
    mp_real_set_fmpz(y, fmpq_numref(a));
    mp_real_set_fmpz(t, fmpq_denref(a));
    mp_real_div(y, y, t, n);
    mp_real_clear(t);
}

/* e^-y y^(e/2) for e = +-1 */
static void
_prefactor(mp_real_t res, const fmpq_t a, int e, slong n)
{
    mp_real_t y, t;
    mp_real_init(y);
    mp_real_init(t);
    _y_ball(y, a, n + 1);
    if (e > 0)
        mp_real_sqrt(t, y, n + 1);
    else
        mp_real_rsqrt(t, y, n + 1);
    mp_real_neg(y, y);
    mp_real_exp_bits(res, y, FLINT_BITS * (n + 1));
    mp_real_mul(res, res, t, n);
    mp_real_clear(y);
    mp_real_clear(t);
}

/* asymptotic series *******************************************************/

/* x as a signed mpn integer; small values are stored in *buf */
static void
_int_set_fmpz(mp_real_hypgeom_int_struct * x, ulong * buf, const fmpz_t c)
{
    if (!COEFF_IS_MPZ(*c))
    {
        *buf = FLINT_UABS(*c);
        x->d = buf;
        x->n = (*buf != 0);
        x->neg = (*c < 0);
    }
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*c);
        x->d = z->_mp_d;
        x->n = FLINT_ABS(z->_mp_size);
        x->neg = (z->_mp_size < 0);
    }
}


typedef struct
{
    const fmpz * p;
    const fmpz * q;
}
_pq_ctx;

/* the number of terms of the asymptotic series for a relative error of
   2^-W (a double estimate; 0 if the smallest term is too large) */
static slong
_asymp_terms(double yd, slong W)
{
    double l = 0.0;
    slong j;

    for (j = 0; ; j++)
    {
        if (l < -(double) W - 8)
            return FLINT_MAX(j, 1);
        if (j + 0.5 >= yd)
            return 0;
        l += log2((j + 0.5) / yd);
    }
}

/* |R_N| <= (1/2)_N y^-N <= N! y^-N (Gamma(N + 1/2) <= N! sqrt(pi)), for
   every N >= 0 (DLMF 8.11(ii)) */
static void
_asymp_bound(mag_t b, const fmpq_t a, slong N)
{
    mag_t t, u;

    mag_init(t);
    mag_init(u);
    mag_fac_ui(b, N);
    mag_set_fmpz_lower(t, fmpq_numref(a));
    mag_set_fmpz(u, fmpq_denref(a));
    mag_div_lower(t, t, u);
    mag_pow_ui_lower(t, t, N);
    mag_div(b, b, t);
    mag_clear(t);
    mag_clear(u);
}

static int
_gamma_asymp(mp_real_t res, const fmpq_t a, slong N, slong n)
{
    mp_real_t s, x0;
    mag_t b;
    slong W = FLINT_BITS * n + 16;
    int ok;

    if (N <= 0)
    {
        N = _asymp_terms(fmpq_get_d(a), W);
        if (N == 0)
            return 0;
    }

    mag_init(b);
    mp_real_init(s);
    mp_real_init(x0);

    /* the remainder after N terms, for every N >= 0 */
    _asymp_bound(b, a, N);
    if (mag_cmp_2exp_si(b, -(W - 8)) > 0)
    {
        /* not accurate enough */
        ok = 0;
        goto cleanup;
    }

    /* sum_(j<N) (-1)^j (1/2)_j y^-j = sum_(k=1)^N prod_(j=1)^(k-1) R(j)/Q(j)
       with Q = 2p, R(j) = -(2j - 1) q = q - 2q j, and P = Q; the sum may
       take a few more terms (N' < 1.09 N), the bound applying to any N */
    {
        mp_real_hypgeom_int_struct P[1], Q[1], R[2];
        ulong buf[3];
        fmpz_t p2, q2, mq2;

        fmpz_init(p2);
        fmpz_init(q2);
        fmpz_init(mq2);
        fmpz_mul_2exp(p2, fmpq_numref(a), 1);
        fmpz_set(q2, fmpq_denref(a));
        fmpz_mul_si(mq2, fmpq_denref(a), -2);
        _int_set_fmpz(P + 0, buf + 0, p2);
        Q[0] = P[0];
        _int_set_fmpz(R + 0, buf + 1, q2);
        _int_set_fmpz(R + 1, buf + 2, mq2);

        N = mp_real_hypgeom_sum(s, x0, P, 1, Q, 1, R, 2, N, 0, n + 2);
        mp_real_div(s, s, x0, n + 2);

        fmpz_clear(p2);
        fmpz_clear(q2);
        fmpz_clear(mq2);
    }

    _asymp_bound(b, a, N);
    if (mag_cmp_2exp_si(b, -(W - 8)) > 0)
    {
        ok = 0;
        goto cleanup;
    }
    _mp_real_add_error_mag(s, b);

    _prefactor(res, a, -1, n + 1);
    mp_real_mul(res, res, s, n);
    ok = 1;

cleanup:
    mp_real_clear(s);
    mp_real_clear(x0);
    mag_clear(b);
    return ok;
}

/* continued fraction *******************************************************/

/* the transposed step matrix [b_k a_k; 1 0], k >= 1 */
static void
_cf_leaf(fmpz_mat_t T, fmpz_t phi, fmpz_t den, slong k, const void * vctx)
{
    const _pq_ctx * c = (const _pq_ctx *) vctx;

    if (k % 2 == 0)
        fmpz_mul_2exp(fmpz_mat_entry(T, 0, 0), c->p, 1);
    else
        fmpz_set_ui(fmpz_mat_entry(T, 0, 0), 2);
    fmpz_mul_ui(fmpz_mat_entry(T, 0, 1), c->q, 2 * (ulong) k);
    fmpz_one(fmpz_mat_entry(T, 1, 0));
    fmpz_one(den);
}

static int
_gamma_cf(mp_real_t res, const fmpq_t a, slong N, slong n)
{
    _mp_real_linrec_struct rec;
    _pq_ctx ctx;
    /* the matrix [R0 R2; R1 R3] (row-major), see below */
    mp_real_struct Rm[4];
    mp_real_ptr r0 = Rm + 0, r2 = Rm + 1, r1 = Rm + 2, r3 = Rm + 3;
    mp_real_t xa, xb, b0, t;
    mag_t m;
    slong W = FLINT_BITS * n + 32, wn = n + 2;
    int ok;

    if (N < 1)
        N = 1;

    ctx.p = fmpq_numref(a);
    ctx.q = fmpq_denref(a);
    rec.d = 2;
    rec.s = 0;
    rec.g = 0;
    rec.leaf = _cf_leaf;
    rec.bottom = NULL;
    rec.bottom_len = 0;
    rec.ctx = &ctx;

    mp_real_init(r0); mp_real_init(r1); mp_real_init(r2); mp_real_init(r3);
    mp_real_init(xa); mp_real_init(xb); mp_real_init(b0); mp_real_init(t);
    mag_init(m);

    /* the transpose of prod_(k=1)^N [b_k 1; a_k 0] =: [R0 R1; R2 R3],
       applied to the identity (all entries are positive) */
    mp_real_set_ui(r0, 1);
    mp_real_zero(r2);
    mp_real_zero(r1);
    mp_real_set_ui(r3, 1);
    _mp_real_linrec_eval(Rm, NULL, NULL, &rec, 1, N + 1, NULL, 2, mp_real_prec_bits(W + 32));

    /* the denominators must be accurate */
    ok = (mp_real_rel_radius_lt_2exp_si(r0) < -64)
        && (mp_real_rel_radius_lt_2exp_si(r1) < -64);
    if (!ok)
        goto cleanup;

    /* the left factor [b0 1; 1 0], b0 = 2p: the last two convergents
       (b0 R0 + R2)/R0 and (b0 R1 + R3)/R1 enclose 2 q x */
    mp_real_set_fmpz(b0, fmpq_numref(a));
    mp_real_mul_2exp_si(b0, b0, 1);
    mp_real_mul(xa, b0, r0, wn);
    mp_real_add(xa, xa, r2, wn);
    mp_real_div(xa, xa, r0, wn);
    mp_real_mul(xb, b0, r1, wn);
    mp_real_add(xb, xb, r3, wn);
    mp_real_div(xb, xb, r1, wn);

    /* the union: the midpoint and half the distance */
    mp_real_sub(t, xa, xb, wn);
    _mp_real_get_mag(m, t);
    mag_mul_2exp_si(m, m, -1);
    mp_real_add(xa, xa, xb, wn);
    mp_real_mul_2exp_si(xa, xa, -1);
    _mp_real_add_error_mag(xa, m);
    if (mp_real_rel_radius_lt_2exp_si(xa) >= -64)
    {
        ok = 0;
        goto cleanup;
    }

    /* Gamma(1/2, y) = e^-y sqrt(y) 2q / (2 q x) */
    _prefactor(res, a, 1, wn);
    mp_real_set_fmpz(t, fmpq_denref(a));
    mp_real_mul_2exp_si(t, t, 1);
    mp_real_mul(res, res, t, wn);
    mp_real_div(res, res, xa, n);

    ok = (mp_real_rel_radius_lt_2exp_si(res) <= -(FLINT_BITS * n - 8));

cleanup:
    mp_real_clear(r0); mp_real_clear(r1); mp_real_clear(r2); mp_real_clear(r3);
    mp_real_clear(xa); mp_real_clear(xb); mp_real_clear(b0); mp_real_clear(t);
    mag_clear(m);
    return ok;
}

/* lower series *************************************************************/

static void
_gamma_lower(mp_real_t res, const fmpq_t a, slong n)
{
    mp_real_hypgeom_series_struct s;
    mp_real_hypgeom_int_struct P[2], R[1];
    ulong buf[5];
    fmpz_t q2, p2;
    mp_real_t sum, t;
    double yd = fmpq_get_d(a);
    slong wn;


    /* the cancellation in sqrt(pi) - gamma(1/2, y) is about e^y */
    wn = n + 3 + (slong) (yd * 1.4426950408889634 / FLINT_BITS);

    fmpz_init(q2);
    fmpz_init(p2);
    mp_real_init(sum);
    mp_real_init(t);

    /* sum_(k>=1) t_k, t_k = prod_(i=1)^(k-1) 2p / (q (2i + 1)):
       P(k) = Q(k) = q (2k + 1), R(k) = 2p */
    fmpz_mul_2exp(q2, fmpq_denref(a), 1);
    fmpz_mul_2exp(p2, fmpq_numref(a), 1);
    _int_set_fmpz(P + 0, buf + 0, fmpq_denref(a));
    _int_set_fmpz(P + 1, buf + 1, q2);
    _int_set_fmpz(R + 0, buf + 2, p2);
    buf[3] = 1;
    buf[4] = 0;
    s.power = 1;
    s.coefP.d = buf + 3; s.coefP.n = 1; s.coefP.neg = 0;
    s.coefD.d = buf + 3; s.coefD.n = 1; s.coefD.neg = 0;
    s.coefQ.d = buf + 4; s.coefQ.n = 0; s.coefQ.neg = 0;
    s.P = P;
    s.Plen = 2;
    s.Q = P;
    s.Qlen = 2;
    s.R = R;
    s.Rlen = 1;
    mp_real_hypgeom_series(sum, &s, wn);

    /* sqrt(pi) - 2 sqrt(y) e^-y sum */
    _prefactor(t, a, 1, wn);
    mp_real_mul(sum, sum, t, wn);
    mp_real_mul_2exp_si(sum, sum, 1);
    mp_real_const_pi4(t, wn, 1);
    mp_real_mul_2exp_si(t, t, 2);
    mp_real_sqrt(t, t, wn);
    /* at the full working precision: the operands are truncated to the
       precision relative to their own magnitude */
    mp_real_sub(res, t, sum, wn);
    mp_real_set_trunc(res, res, n);

    fmpz_clear(q2);
    fmpz_clear(p2);
    mp_real_clear(sum);
    mp_real_clear(t);
}

/* automatic choice *********************************************************/

/* the lower series' term count (an estimate): terms
   (2y)^j / (3 5 ... (2j+1)) relative to their sum, about e^y, below
   2^-(W + 1.44 y) */
static double
_lower_terms(double yd, slong W)
{
    double l = 0.0, target = -(double) W - 1.4426950408889634 * yd - 8;
    slong j;

    /* called only below the asymptotic regime, y < about W, so the loop
       is short (about 3y iterations) */
    for (j = 0; j < WORD_MAX / 4; j++)
    {
        if (j > 2 * yd && l < target)
            return (double) j;
        l += log2(2.0 * yd / (2.0 * j + 3.0));
    }
    return 1e300;
}

static void
_gamma_auto(mp_real_t res, const fmpq_t a, slong n)
{
    double yd;
    slong W = FLINT_BITS * n + 16;

    if (_gamma_asymp(res, a, 0, n))
        return;

    yd = fmpq_get_d(a);

    /* the continued fraction when its estimated length is clearly
       shorter than the lower series' (its matrices cost a bit more) */
    if (yd >= 1.0)
    {
        double ncf = pow(W / 2.4, 2.0) / yd + 1000.0;

        if (1.5 * ncf < _lower_terms(yd, W))
        {
            slong N = (slong) ncf;

            for (;;)
            {
                if (_gamma_cf(res, a, N, n))
                    return;
                if (N > 4 * ncf)
                    break;
                N += N / 2;
            }
        }
    }

    _gamma_lower(res, a, n);
}

/* the public interface: y = p/q from limb arrays *************************/

/* y = p/q in lowest terms; throws unless 0 < y < 2^(FLINT_BITS - 6) (the
   result exponent, about -y log2(e) bits, must stay within the range of
   mp_real_exp) */
static void
_y_set(fmpq_t y, nn_srcptr p, slong pn, nn_srcptr q, slong qn)
{
    fmpz_t t;
    int big;

    while (pn > 0 && p[pn - 1] == 0)
        pn--;
    while (qn > 0 && q[qn - 1] == 0)
        qn--;
    if (pn <= 0 || qn <= 0)
        flint_throw(FLINT_ERROR, "mp_real_gamma_upper_half_frac: p = 0 or q = 0\n");

    fmpz_set_ui_array(fmpq_numref(y), p, pn);
    fmpz_set_ui_array(fmpq_denref(y), q, qn);
    fmpq_canonicalise(y);

    fmpz_init(t);
    fmpz_mul_2exp(t, fmpq_denref(y), FLINT_BITS - 6);
    big = fmpz_cmp(fmpq_numref(y), t) >= 0;
    fmpz_clear(t);
    if (big)
        flint_throw(FLINT_ERROR, "mp_real_gamma_upper_half_frac: "
            "y >= 2^%wd\n", (slong) (FLINT_BITS - 6));
}

void
mp_real_gamma_upper_half_frac(mp_real_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, slong n)
{
    fmpq_t y;
    fmpq_init(y);
    _y_set(y, p, pn, q, qn);
    _gamma_auto(res, y, n);
    fmpq_clear(y);
}

int
_mp_real_gamma_upper_half_frac_asymp(mp_real_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, slong N, slong n)
{
    fmpq_t y;
    int ok;
    fmpq_init(y);
    _y_set(y, p, pn, q, qn);
    ok = _gamma_asymp(res, y, N, n);
    fmpq_clear(y);
    return ok;
}

int
_mp_real_gamma_upper_half_frac_cf(mp_real_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, slong N, slong n)
{
    fmpq_t y;
    int ok;
    fmpq_init(y);
    _y_set(y, p, pn, q, qn);
    ok = _gamma_cf(res, y, N, n);
    fmpq_clear(y);
    return ok;
}

void
_mp_real_gamma_upper_half_frac_lower(mp_real_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, slong n)
{
    fmpq_t y;
    fmpq_init(y);
    _y_set(y, p, pn, q, qn);
    _gamma_lower(res, y, n);
    fmpq_clear(y);
}
