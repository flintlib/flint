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
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* the documented maximum of *err, checked on export */
#define NEGLOG_NEWTON_MAX_ERR 2
#define ATAN_NEWTON_MAX_ERR 2

/* Newton-Taylor evaluation of -log x on [1/2, 1) and atan x on [0, 1)
   from the forward functions (fixed-point ports of arb_log_newton and
   arb_atan_newton).

   THE STEP.  Let f be -log or atan, p = FLINT_BITS wn the working
   precision and t an approximation of f(x) correct to about p0 bits.
   The forward function evaluated at t at full precision gives, in
   exact arithmetic, a small w with f(x) = t + g(w) for g the Taylor
   series of log(1 + w) or atan(w):

     -log x:  E = exp(t) in [1, 2), w = x E - 1,
              -log x = t - log(1 + w)                (one product);
     atan x:  (s, c) = (sin t, cos t),
              w = (x c - s) / (x s + c),
              atan x = t + atan(w)                    (two products
              and a division at (p - p0) bits, the top p0 bits of the
              quotient being zero).

   These identities hold exactly for the fixed-point number t, whatever
   its error: t only has to be good enough for w to be small, and the
   inaccuracy of the forward function enters the result only through
   its (small) effect on w.  With |w| < 2^-p0, d + 1 terms of g suffice
   when (d + 1) p0 >= p, and since the k-th term has k p0 leading zero
   bits the polynomial costs a fraction of a multiplication at
   precision p.  The whole step is therefore one forward evaluation
   plus a recursive f at 1/(d + 1) of the precision: below the
   MP_REAL_NEWTON_CUTOFF the starting value comes from the bitwise
   functions, above it from this step itself, so the recursion costs a
   geometric series worth about 1/d of a forward evaluation.

   THE ARITHMETIC is mp_real: the forward function (a fixed-point
   routine, imported with its documented error) and the starting
   value (the bitwise inverse, or this step itself with its midpoint
   truncated to p0 bits -- the identities hold for any fixed-point t)
   are the only mpn inputs.  The residual w = x E - 1 resp.
   (x c - s) / (x s + c) is ball arithmetic, with the cancellation
   and the forward errors carried by the balls (the atan denominator
   and the division at the precision the cancelled numerator has
   left); the series g(w) = w P(v), v = w (log) or w^2 (atan), has
   P(v) = sum_{k<N} c_k v^k / den with the Taylor coefficients c_k
   scaled to integers by the lcm den of their denominators (below
   2^31), its powers of v by squaring at the precision each term
   contributes at (the k-th has k p0 leading zero bits), one
   mp_real_mul_ui per term, one product by w and one mp_real_div_ui by
   den; the tail of the series is added to the radius from the
   rigorous bound |w| < 2^uexp of the ball (and is below a quarter of
   a working ulp by the choice of p0, see _newton_start_limbs), and
   the result is t -+ correction as a ball.  The export to n limbs
   truncates (one ulp) and checks the ball's radius -- about 2^-40
   ulps: the forward errors (3 or 4 working ulps of B^-(n+1)) enter
   through w only -- against NEGLOG_NEWTON_MAX_ERR =
   ATAN_NEWTON_MAX_ERR = 2.

   An earlier implementation in hand-written mpn arithmetic (windowed
   middle products of the top limbs, a rectangular splitting of the
   polynomial with mpn_addmul_1 rows, and an error accounting of its
   own) measured the same speed at every size from 100 to 64000
   limbs -- the forward function is 90% of the time -- and was
   dropped.

   The sums t +- corr never underflow: -log x >= 1 - x >= B^-n is a
   full B working ulps, far above the errors, and atan adds t >= 0
   and a correction whose sign it knows. */

/* extra bits of margin between the starting value's error and 2^-p0
   in the truncation bound: the bitwise functions are within 2^11.1
   ulps of their precision, x differs from its top p0n limbs by 2 ulps
   at most, the recursive start is within 2 ulps */
#define NEWTON_SLACK_BITS 16

/* precision (limbs) of the starting value for a series truncated
   after degree d: (d + 1)(FLINT_BITS p0n - NEWTON_SLACK_BITS) >=
   FLINT_BITS wn + 2 leaves the tail |w|^(d+1) / (d + 1) below a
   quarter of a working ulp */
/* the top p0n limbs of (x, n), zero-extended below when p0n > n (on
   32-bit machines the bitwise kernels of the starting values need two
   limbs, and p0n = 2 > n for n = 1); buf holds p0n limbs */
static nn_srcptr
_newton_start_top(nn_ptr buf, nn_srcptr x, slong n, slong p0n)
{
    if (p0n <= n)
        return x + (n - p0n);
    flint_mpn_zero(buf, p0n - n);
    flint_mpn_copyi(buf + (p0n - n), x, n);
    return buf;
}

#if FLINT_BITS == 64
#define NEWTON_START_MIN 1
#else
#define NEWTON_START_MIN 2
#endif

static slong
_newton_start_limbs(slong wn, slong d)
{
    slong num = FLINT_BITS * wn + 2 + (d + 1) * NEWTON_SLACK_BITS;
    return (num + FLINT_BITS * (d + 1) - 1) / (FLINT_BITS * (d + 1));
}

/* Taylor coefficients: log(1 + w) = w sum_{k>=0} (-1)^k w^k / (k + 1)
   and atan(w) = w sum_{k>=0} (-1)^k w^(2k) / (2k + 1), as magnitudes
   scaled by the lcm of the denominators up to the maximal order:
   lcm(1, ..., 16) = 720720 and lcm(1, 3, ..., 25) = 1673196525, both
   below 2^31 so that they are words on 32-bit limbs too (one more
   atan term would need the factors 27, 29, 31 and a 43-bit lcm) */
#define NEWTON_MAX_N 16
#define NEWTON_LOG_MAX_N 16
#define NEWTON_ATAN_MAX_N 13
static const ulong newton_log_coeffs[NEWTON_LOG_MAX_N] = {
    720720, 360360, 240240, 180180, 144144, 120120, 102960, 90090,
    80080, 72072, 65520, 60060, 55440, 51480, 48048, 45045 };
#define NEWTON_LOG_DEN UWORD(720720)
static const ulong newton_atan_coeffs[NEWTON_ATAN_MAX_N] = {
    UWORD(1673196525), UWORD(557732175), UWORD(334639305), UWORD(239028075),
    UWORD(185910725), UWORD(152108775), UWORD(128707425), UWORD(111546435),
    UWORD(98423325), UWORD(88062975), UWORD(79676025), UWORD(72747675),
    UWORD(66927861) };
#define NEWTON_ATAN_DEN UWORD(1673196525)

/* the default numbers of terms (measured, see tune/tune-newton.c):
   the total is flat to within a few percent over a wide range
   around these, more terms trading a cheaper starting value (at
   1/(N + 1) resp. 1/(2N) of the precision) for a longer polynomial,
   which pays off slowly as the precision grows */
static slong
_neglog_default_N(slong n)
{
    return (n <= 64) ? 6 : (n <= 2000) ? 8 : (n <= 8000) ? 10 : 12;
}

static slong
_atan_default_N(slong n)
{
    return (n <= 64) ? 3 : (n <= 2000) ? 4 : (n <= 8000) ? 5 : 6;
}

/* ==== the steps ==========================================================*/

/* res = sum_{k=1}^{N} s_k (den / d_k) v^(k-1) w for the series
   w sum_k s_k v^(k-1) / d_k with v = w (log: d_k = k) or v = w^2
   (atan: d_k = 2k - 1), s_k = (-1)^(k+1) when alternating; the
   powers of v by squaring at the precision they contribute at, from
   uexp with |w| < 2^uexp; W accumulated at the precision of the
   leading term; the caller divides by den and adds the tail */
static void
_mp_real_newton_series(mp_real_t res, const mp_real_t w, slong uexp, slong p,
    const ulong * c, slong N, int odd, int alternating)
{
    mp_real_struct pw[NEWTON_MAX_N + 1];
    mp_real_t v, W, T;
    slong k, cp, vexp = odd ? 2 * uexp : uexp;
    slong cw = FLINT_MAX(2, p + uexp / FLINT_BITS + 2);   /* W is then multiplied by w */

    mp_real_init(v);
    mp_real_init(W);
    mp_real_init(T);

    /* the sum over v: W = sum_{k=1}^N s_k c_k v^(k-1) */
    mp_real_set_ui(W, c[0]);
    if (N > 1)
    {
        if (odd)
            mp_real_mul(v, w, w, FLINT_MAX(2, p + vexp / FLINT_BITS + 2));
        else
            mp_real_set(v, w);
        for (k = 1; k < N; k++)
        {
            const mp_real_struct * vk;

            /* v^k needed to p + (k vexp + uexp) / B limbs */
            cp = FLINT_MAX(2, p + (k * vexp + uexp) / FLINT_BITS + 2);
            if (k == 1)
                vk = v;
            else
            {
                mp_real_init(pw + k);
                if (k % 2 == 0)
                    mp_real_mul(pw + k, (k / 2 == 1) ? v : pw + k / 2,
                        (k / 2 == 1) ? v : pw + k / 2, cp);
                else
                    mp_real_mul(pw + k, pw + k - 1, v, cp);
                vk = pw + k;
            }
            if (alternating && (k & 1))
                mp_real_submul_ui(W, W, vk, c[k], cw);
            else
                mp_real_addmul_ui(W, W, vk, c[k], cw);
        }
        for (k = 2; k < N; k++)
            mp_real_clear(pw + k);
    }

    mp_real_mul(res, W, w, p);

    mp_real_clear(v);
    mp_real_clear(W);
    mp_real_clear(T);
}

/* -log(x) as a ball for (x, n) in [1/2, 1); forward, N as for the
   tunable fixed-point worker */
void
_mp_real_neglog_newton_ball(mp_real_t res, nn_srcptr x, slong n, int forward,
    slong N)
{
    slong wn, p0n, uexp;
    ulong ferr;
    nn_ptr t, E;
    nn_srcptr xtop;
    mp_real_t xb, tb, w, corr;
    TMP_INIT;

    FLINT_ASSERT(n >= 1 && (x[n - 1] >> (FLINT_BITS - 1)) != 0);

    if (N == 0)
        N = _neglog_default_N(n);
    N = FLINT_MIN(N, NEWTON_LOG_MAX_N);

    wn = n + 1;
    p0n = FLINT_MIN(_newton_start_limbs(wn, N), n);
    p0n = FLINT_MAX(p0n, NEWTON_START_MIN);     /* <= wn */

    TMP_START;
    t = TMP_ALLOC(wn * sizeof(ulong));
    /* wn + 1 limbs below, and the starting value's log 2 and argument
       side by side (p0n limbs each, p0n <= wn) */
    E = TMP_ALLOC(FLINT_MAX(wn + 1, 2 * p0n) * sizeof(ulong));
    xtop = _newton_start_top(TMP_ALLOC(p0n * sizeof(ulong)), x, n, p0n);
    mp_real_init(xb);
    mp_real_init(tb);
    mp_real_init(w);
    mp_real_init(corr);

    /* the starting value t = -log(x_top) at p0n limbs, exact as a
       fixed-point number */
    flint_mpn_zero(t, wn - p0n);
    if (p0n <= MP_REAL_NEWTON_CUTOFF)
    {
        nn_ptr l2 = E, arg = E + p0n;
        _mp_real_const_log2(l2, NULL, p0n, 1);
        mpn_lshift(arg, xtop, p0n, 1);   /* 2 x_top - 1 */
        _mp_real_log1p_bitwise_rs(t + (wn - p0n), NULL, arg, p0n, 0);
        if (mpn_sub_n(t + (wn - p0n), l2, t + (wn - p0n), p0n))
            flint_mpn_zero(t + (wn - p0n), p0n);
    }
    else
    {
        _mp_real_neglog_newton_ball(corr, xtop, p0n, forward, 0);
        _mp_real_get_fixed(t + (wn - p0n), NULL, corr, p0n);
    }
    _mp_real_set_mpn_2exp(tb, t, wn, -FLINT_BITS * wn);

    /* E = exp(t) in [1, 2), imported with its error */
    if (forward == 2)
        _mp_real_exp_notab(E, &ferr, t, wn);
    else if (forward)
        _mp_real_exp_bitwise_rs(E, &ferr, t, wn, 0);
    else
        _mp_real_exp_diophantine(E, &ferr, t, wn);
    _mp_real_set_mpn_2exp(w, E, wn + 1, -FLINT_BITS * wn);
    _mp_real_add_error_ulps_at(w, ferr, -wn);

    /* w = x E - 1 */
    _mp_real_set_mpn_2exp(xb, x, n, -FLINT_BITS * n);
    mp_real_mul(w, w, xb, wn + 1);
    mp_real_set_ui(xb, 1);
    mp_real_sub(w, w, xb, wn + 1);
    uexp = mp_real_abs_bound_lt_2exp_si(w);

    /* -log x = t - log(1 + w), the series to N terms, its tail
       |w|^(N+1) / ((N + 1)(1 - |w|)) below 2^((N+1) uexp + 1) */
    if (uexp > -2)
        flint_throw(FLINT_ERROR, "_mp_real_neglog_newton: residual too large\n");
    _mp_real_newton_series(corr, w, uexp, wn + 1, newton_log_coeffs, N, 0, 1);
    mp_real_div_ui(corr, corr, NEWTON_LOG_DEN, wn + 1);
    mp_real_add_error_2exp_si(corr, (N + 1) * uexp + 1);
    mp_real_sub(res, tb, corr, wn + 1);

    mp_real_clear(xb);
    mp_real_clear(tb);
    mp_real_clear(w);
    mp_real_clear(corr);
    TMP_END;
}

/* atan(x) as a ball for (x, n) in [0, 1) */
void
_mp_real_atan_newton_ball(mp_real_t res, nn_srcptr x, slong n, int forward,
    slong N)
{
    slong wn, p0n, uexp, wq;
    ulong ferr;
    nn_ptr t, s, c;
    nn_srcptr xtop;
    mp_real_t xb, tb, sb, cb, num, den, w, corr;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);

    if (N == 0)
        N = _atan_default_N(n);
    N = FLINT_MIN(N, NEWTON_ATAN_MAX_N);

    wn = n + 1;
    p0n = FLINT_MIN(_newton_start_limbs(wn, 2 * N - 1), n);
    p0n = FLINT_MAX(p0n, NEWTON_START_MIN);     /* <= wn */

    TMP_START;
    t = TMP_ALLOC(wn * sizeof(ulong));
    xtop = _newton_start_top(TMP_ALLOC(p0n * sizeof(ulong)), x, n, p0n);
    s = TMP_ALLOC((wn + 1) * sizeof(ulong));
    c = TMP_ALLOC((wn + 1) * sizeof(ulong));
    mp_real_init(xb);
    mp_real_init(tb);
    mp_real_init(sb);
    mp_real_init(cb);
    mp_real_init(num);
    mp_real_init(den);
    mp_real_init(w);
    mp_real_init(corr);

    /* the starting value t = atan(x_top) at p0n limbs */
    flint_mpn_zero(t, wn - p0n);
    if (p0n <= MP_REAL_NEWTON_CUTOFF)
        _mp_real_atan_bitwise_rs(t + (wn - p0n), NULL, xtop, p0n, 0);
    else
    {
        _mp_real_atan_newton_ball(corr, xtop, p0n, forward, 0);
        _mp_real_get_fixed(t + (wn - p0n), NULL, corr, p0n);
    }
    _mp_real_set_mpn_2exp(tb, t, wn, -FLINT_BITS * wn);

    /* (s, c) = (sin t, cos t), imported with their errors */
    if (forward == 2)
        _mp_real_sin_cos_notab(s, c, &ferr, t, wn);
    else if (forward)
        _mp_real_sin_cos_bitwise_rs(s, c, &ferr, t, wn, 0);
    else
        _mp_real_sin_cos_diophantine(s, c, &ferr, t, wn);
    _mp_real_set_mpn_2exp(sb, s, wn + 1, -FLINT_BITS * wn);
    _mp_real_set_mpn_2exp(cb, c, wn + 1, -FLINT_BITS * wn);
    _mp_real_add_error_ulps_at(sb, ferr, -wn);
    _mp_real_add_error_ulps_at(cb, ferr, -wn);

    /* w = (x c - s) / (x s + c): the numerator cancels to about
       wn - p0n limbs, and the denominator (>= 1) is only needed to
       that precision */
    _mp_real_set_mpn_2exp(xb, x, n, -FLINT_BITS * n);
    mp_real_mul(num, xb, cb, wn + 1);
    mp_real_sub(num, num, sb, wn + 1);
    wq = FLINT_MAX(2, wn - p0n + 3);
    mp_real_mul(den, xb, sb, wq);
    mp_real_add(den, den, cb, wq);
    mp_real_div(w, num, den, wq);
    uexp = mp_real_abs_bound_lt_2exp_si(w);

    /* atan x = t + atan(w), the series in w^2 to N terms, its tail
       |w|^(2N+1) / ((2N + 1)(1 - w^2)) below 2^((2N+1) uexp + 1) */
    if (uexp > -2)
        flint_throw(FLINT_ERROR, "_mp_real_atan_newton: residual too large\n");
    _mp_real_newton_series(corr, w, uexp, wn + 1, newton_atan_coeffs, N, 1, 1);
    mp_real_div_ui(corr, corr, NEWTON_ATAN_DEN, wn + 1);
    mp_real_add_error_2exp_si(corr, (2 * N + 1) * uexp + 1);
    mp_real_add(res, tb, corr, wn + 1);

    mp_real_clear(xb);
    mp_real_clear(tb);
    mp_real_clear(sb);
    mp_real_clear(cb);
    mp_real_clear(num);
    mp_real_clear(den);
    mp_real_clear(w);
    mp_real_clear(corr);
    TMP_END;
}

/* the fixed-point exports of the ball steps, checked against the
   budgets; forward 0 = the diophantine forward function, 1 = the
   bitwise one, 2 = the table-free one (exp_notab, sin_cos_notab), N = 0 the default number of terms */
void
_mp_real_neglog_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward,
    slong N)
{
    mp_real_t r;
    ulong bound;
    mp_real_init(r);
    _mp_real_neglog_newton_ball(r, x, n, forward, N);
    _mp_real_get_fixed(y, &bound, r, n);
    if (!(bound <= NEGLOG_NEWTON_MAX_ERR))
        flint_throw(FLINT_ERROR, "_mp_real_neglog_newton: error bound %wu ulps\n", bound);
    if (err != NULL)
        *err = bound;
    mp_real_clear(r);
}

void
_mp_real_atan_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward,
    slong N)
{
    mp_real_t r;
    ulong bound;
    mp_real_init(r);
    _mp_real_atan_newton_ball(r, x, n, forward, N);
    _mp_real_get_fixed(y, &bound, r, n);
    if (!(bound <= ATAN_NEWTON_MAX_ERR))
        flint_throw(FLINT_ERROR, "_mp_real_atan_newton: error bound %wu ulps\n", bound);
    if (err != NULL)
        *err = bound;
    mp_real_clear(r);
}

void
_mp_real_neglog_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_neglog_newton_tune(y, err, x, n, 0, 0);
}

void
_mp_real_atan_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_atan_newton_tune(y, err, x, n, 0, 0);
}
