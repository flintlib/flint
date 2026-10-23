/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* the documented maximum of *err, checked on export */
#define NEGLOG_NEWTON_MAX_ERR 2

/* -log(x) for x in [1/2, 1) by the Sasaki-Kanada formula

       log(1/q) = pi / agm(theta_2(q)^2, theta_3(q)^2),   0 < q < 1,

   which is exact: only the theta series are truncated.  With
   r = x 2^-e and q = r^4 (so that theta_2(q) = 2 r S_2 needs no
   fourth root),

       theta_3 = 1 + 2 sum_{n>=1} q^(n^2),
       S_2 = 1 + sum_{n>=1} q^(n(n+1)),
       -log x = (pi/4) / agm(4 r^2 S_2^2, theta_3^2) - e log 2,

   whose first AGM step needs no square root (see below).

   THE THETA SERIES.  The exponents of the two series merged in order
   are floor(m^2/4), m = 2, 3, ... (1, 2, 4, 6, 9, 12, 16, ...: even m
   for theta_3, odd m for S_2), and consecutive ones differ by
   floor(m/2): each term is the previous one times q^j, j = floor(m/2),
   and the powers q^j come by squaring for even j and one product by q
   for odd j, each to the precision of the term q^(j^2) that first uses
   it.  So the two series cost one product per term and one squaring
   or product per two terms, each at the precision the term
   contributes at (q^k is 2^(-4ek) relative to the leading 1).  The terms kept are those with
   4 e k below the working precision; the omitted ones, of exponents
   k >= k_next in each series and ratios below q <= 1/2, sum to less
   than 2 q^k_next, added to the radius.

   THE BALANCE.  e = ceil(p / (4 (N + 1)^2)) for N terms of theta_3
   (and N of S_2), p the working precision in bits: the AGM of
   (4 r^2 S_2^2 ~ 2^(2 - 2e), ~1) takes about log2(p) + log2(e)
   iterations -- about log2(e) of them in the slow phase before the
   quadratic convergence sets in -- so each doubling of N saves two
   iterations (a product and a square root each) for about 3N more
   products at decreasing precision; the default N is tuned below.

   ACCURACY.  The two terms are about e log 2 in size and the result
   below log 2, so the working precision carries one guard limb for
   the cancellation (e < 2^FLINT_BITS); pi/4 and log 2 come from the
   per-thread caches of verified floors (within one ulp), and the ball
   arithmetic carries all errors to a rigorous bound, checked on
   export. */

/* the default number of theta_3 terms: measured flat to within a
   few percent from N = 1 to 6 at every size from 64 to 65536 limbs
   (the merged series cost about 3N products at an average of 0.6 of
   the precision, and each doubling of N saves two AGM iterations of a
   product and a square root each), slower from 8 on */
static slong
_log_agm_default_N(slong n)
{
    return 4;
}

void
_mp_real_neglog_agm_ball(mp_real_t res, nn_srcptr x, slong n, slong N)
{
    slong wp = n + 3, p, e, m, j, k, kn, qe, cp;
    mp_real_t r2, q, t, T3, S2, a, b, c;
    nn_ptr tmp;

    FLINT_ASSERT(n >= 1 && (x[n - 1] >> (FLINT_BITS - 1)) != 0);

    if (N == 0)
        N = _log_agm_default_N(n);
    N = FLINT_MAX(N, 1);
    p = FLINT_BITS * wp;
    e = (p + 4 * (N + 1) * (N + 1) - 1) / (4 * (N + 1) * (N + 1));
    e = FLINT_MAX(e, 2);

    mp_real_init(r2); mp_real_init(q); mp_real_init(t);
    mp_real_init(T3); mp_real_init(S2); mp_real_init(a); mp_real_init(b);
    mp_real_init(c);

    /* r^2 = x^2 2^-2e to full precision, q = r^4 to the precision it
       contributes at (2^(-4e) below the leading 1) */
    _mp_real_set_mpn_2exp(t, x, n, -FLINT_BITS * n);
    mp_real_mul(r2, t, t, wp);
    mp_real_mul_2exp_si(r2, r2, -2 * e);
    mp_real_mul(q, r2, r2, FLINT_MAX(2, wp - (4 * e) / FLINT_BITS + 1));
    qe = mp_real_abs_bound_lt_2exp_si(q);

    /* the merged series: t = q^floor(m^2/4) = t_(m-1) q^floor(m/2),
       with the powers q^j in a table by squaring for even j
       (q^j = (q^(j/2))^2) and one product by q for odd j, each to the
       precision of the term q^(j^2) that first uses it (it only ever
       multiplies terms at least that small); the terms q^2 = q^(2/2)^2
       and q^4 = (q^2)^2 are the power itself and a squaring */
    {
        slong jmax = N + 9;
        mp_real_struct * pw = flint_malloc((jmax + 1) * sizeof(mp_real_struct));

        for (j = 0; j <= jmax; j++)
            mp_real_init(pw + j);
        mp_real_set(pw + 1, q);

        mp_real_zero(T3);
        mp_real_zero(S2);
        mp_real_set(t, q);        /* m = 2, k = 1 */
        j = 1;
        k = 1;
        for (m = 2; ; m++)
        {
            if (m > 2)
            {
                if (m % 2 == 0 && m > 4)
                {
                    /* the power q^(m/2) (q^2 is the term of m = 3) */
                    j = m / 2;
                    cp = FLINT_MAX(2, wp + (j * j * qe) / FLINT_BITS + 1);
                    if (j % 2 == 0)
                        mp_real_mul(pw + j, pw + j / 2, pw + j / 2, cp);
                    else
                        mp_real_mul(pw + j, pw + j - 1, q, cp);
                }
                k += m / 2;
                cp = FLINT_MAX(2, wp + (k * qe) / FLINT_BITS + 1);
                if (m == 3)
                {
                    /* t = q^2 by a squaring, at the precision of the
                       term (k = 2), which also serves as the power */
                    mp_real_mul(pw + 2, q, q, cp);
                    mp_real_set(t, pw + 2);
                }
                else if (m == 4)
                    mp_real_mul(t, pw + 2, pw + 2, cp);       /* q^4 */
                else
                    mp_real_mul(t, t, pw + m / 2, cp);
            }
            if (m % 2 == 0)
                mp_real_add(T3, T3, t, wp);
            else
                mp_real_add(S2, S2, t, wp);

            /* the next exponent; stop when it is below the precision */
            kn = k + (m + 1) / 2;
            if (kn * qe < -p - 8 || m / 2 >= jmax - 1)
                break;
        }

        for (j = 0; j <= jmax; j++)
            mp_real_clear(pw + j);
        flint_free(pw);
    }
    /* the tails of both series: below 2 q^kn */
    mp_real_add_error_2exp_si(T3, kn * qe + 1);
    mp_real_add_error_2exp_si(S2, kn * qe + 1);

    /* the first AGM step without a square root: with u = 2 r S_2 and
       theta_3, the pair (u^2, theta_3^2) goes to
       ((u^2 + theta_3^2)/2, u theta_3) -- the doubling formulas
       theta_3(q)^2 + theta_2(q)^2 = theta_3(q^(1/2))^2 and
       2 theta_2(q) theta_3(q) = theta_2(q^(1/2))^2: the theta
       functions at the square root of the nome, halved */
    mp_real_set_ui(t, 1);
    mp_real_add(S2, S2, t, wp);
    _mp_real_set_mpn_2exp(a, x, n, -FLINT_BITS * n);
    mp_real_mul(S2, S2, a, wp);
    mp_real_mul_2exp_si(S2, S2, 1 - e);           /* u = 2 r S_2 */
    mp_real_mul_2exp_si(T3, T3, 1);
    mp_real_add(T3, T3, t, wp);               /* theta_3 */
    mp_real_mul(b, S2, T3, wp);               /* u theta_3 */
    mp_real_mul(a, S2, S2, wp);
    mp_real_mul(t, T3, T3, wp);
    mp_real_add(a, a, t, wp);
    mp_real_mul_2exp_si(a, a, -1);               /* (u^2 + theta_3^2)/2 */

    mp_real_agm(c, a, b, wp);

    /* (pi/4) / agm - e log 2 */
    tmp = flint_malloc(wp * sizeof(ulong));
    _mp_real_const_pi4(tmp, NULL, wp, 1);
    _mp_real_set_mpn_2exp(t, tmp, wp, -FLINT_BITS * wp);
    _mp_real_add_error_ulps_at(t, 1.0, -wp);
    mp_real_div(c, t, c, wp);
    _mp_real_const_log2(tmp, NULL, wp, 1);
    _mp_real_set_mpn_2exp(t, tmp, wp, -FLINT_BITS * wp);
    _mp_real_add_error_ulps_at(t, 1.0, -wp);
    mp_real_submul_ui(res, c, t, (ulong) e, wp);
    flint_free(tmp);

    mp_real_clear(r2); mp_real_clear(q); mp_real_clear(t);
    mp_real_clear(T3); mp_real_clear(S2); mp_real_clear(a); mp_real_clear(b);
    mp_real_clear(c);
}

void
_mp_real_neglog_agm_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong N)
{
    mp_real_t r;
    ulong bound;
    mp_real_init(r);
    _mp_real_neglog_agm_ball(r, x, n, N);
    _mp_real_get_fixed(y, &bound, r, n);
    if (!(bound <= NEGLOG_NEWTON_MAX_ERR))
        flint_throw(FLINT_ERROR, "_mp_real_neglog_agm: error bound %wu ulps\n", bound);
    if (err != NULL)
        *err = bound;
    mp_real_clear(r);
}

void
_mp_real_neglog_agm(nn_ptr y, ulong * err, nn_srcptr x, slong n)
{
    _mp_real_neglog_agm_tune(y, err, x, n, 0);
}
