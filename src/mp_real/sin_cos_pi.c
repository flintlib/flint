/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "mp_real.h"
#include "impl.h"

/* sin(pi p/q) and cos(pi p/q) for words p, q.

   Closed forms for q <= 6 and (square roots, beating both algorithms
   below from a few limbs) q = 5, 8, 10, 12; otherwise two algorithms:

   KERNEL.  v = pi p/q (or its complement pi/2 - v, so that v <= pi/4)
   formed in a fixed-point buffer from the cached pi/4 (one mpn_mul_1
   and one mpn_divrem_1, two ulps below the exact value), and the
   kernel of mp_real_sin_cos_bits.  Its cost, relative to one
   multiplication M(n), is about 16-25 up to a few hundred limbs and
   grows slowly beyond (about 40 at 1024 limbs, 60 at 4096, 85 at
   16384, 150 at 65536).

   CHEBYSHEV.  x* = cos phi*, phi* = pi R/Q in (0, pi/2), is the simple
   root of the Chebyshev polynomial U_(Q-1) with T_Q(x*) = eps =
   (-1)^R.  Given an exact approximation x0 = cos phi0 and writing
   phi0 = phi* + delta, the values

       C = T_Q(x0) = eps cos(Q delta),
       A = x0 T_Q(x0) - T_(Q+1)(x0) = eps s0 sin(Q delta),  s0 = sin phi0,

   come from the order-doubling (Montgomery ladder) evaluation of the
   pair (T_k, T_(k+1)) -- one squaring and one product per bit of Q, in
   the form V_k = 2 T_k: V_2k = V_k^2 - 2, V_(2k+1) = V_k V_(k+1) - V_1
   (for odd Q = 2k + 1 the last doubling and the product x0 T_Q become
   V_Q = V_k V_(k+1) - V_1 and A = (V_k - V_(k+1))(V_k + V_(k+1))/4).
   With a = eps A, sigma = sin(Q delta) = a / s0, y = sigma^2 =
   (1 - eps C)(1 + eps C) and b = 1/Q, exactly

       x* = cos(phi0 - delta) = x0 cos delta + s0 sin delta
          = x0 + x0 G(y) + (a/Q) F(y),

       G(y) = cos(b asin sigma) - 1 = 2F1(-b/2, b/2; 1/2; y) - 1
            = sum_{j>=1} c_j y^j,
       F(y) = Q sin(b asin sigma) / sigma = 2F1((1-b)/2, (1+b)/2; 3/2; y)
            = sum_{j>=0} f_j y^j,

   with c_1 = -1/(2Q^2), c_(j+1) = c_j (4j^2 Q^2 - 1) / (2Q^2 (j+1)(2j+1))
   and f_0 = 1, f_(j+1) = f_j ((2j+1)^2 Q^2 - 1) / (2Q^2 (2j+3)(j+1)):
   the c_j (j >= 1) are negative and decrease in magnitude from 1/(2Q^2),
   the f_j positive and decreasing from 1, so for |y| <= 1/2 the tails
   after y^J are below |y|^J / Q^2 resp. 2 |y|^J.  No division and no
   square root enters: s0 appears only through a = s0 sigma and y.
   This is the identity when |Q delta| < pi/2 (asin's principal
   branch), which the step checks from a rigorous bound on |x0 - x*|.
   The series share the powers of z = y/(4Q^2), in which both have
   integer coefficients over the common denominators m! (2m-1)!! resp.
   m! (2m+1)!! (m terms): a term costs one product by a short exact
   integer, and each series one mp_real_div_ui.

   With x0 accurate to e bits, sigma ~ Q 2^-e, and summing the series
   until the terms fall below the target gives a step of any order: x0
   at a fraction 1/r of the precision (plus (r-1) log2 Q bits, the
   growth of the residual) and the evaluation at full precision.  The
   ladder dominates, about (log2 Q - log2 r) (S + M)(n), the first
   log2 r doublings of the short x0 being short products; the series
   adds a few products at falling precision (about 3 M(n) at order 12).
   The recursion runs in ball arithmetic (it handles the propagation of
   rounding errors through the ladder, whose radius grows by 2 bits per
   doubling), the radius of each level discarded except as the rigorous
   bound that the next step's branch check reads; it bottoms out in the
   kernel as soon as that is cheaper at the level's precision, which
   makes "kernel at n/r, then one step" the choice at large Q.  The
   other value, sin phi* = sqrt(1 - x*^2), costs one square root, about
   3 M(n) (the short x0 does not make the square root of 1 - x0^2
   cheaper).  Of the two targets cos(pi t/q) and
   sin(pi t/q) = cos(pi (q - 2t)/(2q)), the one with the smaller
   denominator is iterated (its conditioning, sin phi* >= 2R/Q, is paid
   in guard bits), the other value taken by the square root.

   Order r: 12 (measured from 192 to 65536 limbs: 8, 12 and 16 within a
   few percent of each other, 6 about 9% and 4 about 17% slower).

   The complex alternative, the reciprocal Q-th root w = exp(-i phi*)
   of eps by Newton-Taylor steps (w0 (1 - u)^(-1/Q), u = 1 - eps w0^Q,
   as for the real roots), needs a complex squaring (two products) per
   bit and a complex product by the short w0 per set bit where the
   ladder has one squaring and one product; in the same ball
   arithmetic (tune/tune-sin-cos-pi.c) it measured 1.4-2 times slower
   than the Chebyshev iteration with its square root, which makes up
   less than the complex iteration's extra cost even at Q = 5.

   Choice: the Chebyshev iteration for denominators up to a measured
   bit length depending on n (none below 48 limbs; 3 bits at 64 limbs,
   9 at 256, 17 at 1024, 25 at 4096, 40 at 16384, 48 at 65536, all
   admissible Q beyond, where the kernel no longer uses tables). */

/* the largest bit length of the denominator for which the Chebyshev
   iteration beats the kernel at n limbs (tune/tune-sin-cos-pi.c;
   measured with both outputs on 64-bit x86, for which one output moves
   the crossovers by a bit or two; beyond 65536 limbs the kernel
   switches to the untabulated one and the iteration wins for every
   admissible denominator) */
static slong
_scp_cheb_kbmax(slong n)
{
    /* by 64-bit limbs of precision (on 32-bit machines, untuned) */
    static const int tab[][2] = {
        {48, 0}, {96, 3}, {160, 7}, {224, 8}, {320, 9}, {448, 10},
        {640, 12}, {1280, 17}, {1792, 19}, {2560, 20}, {6144, 25},
        {12288, 30}, {24576, 40}, {49152, 43}, {98304, 48}};
    slong i;

    n = n / (64 / FLINT_BITS);
    for (i = 0; i < 15; i++)
        if (n < tab[i][0])
            return tab[i][1];
    return FLINT_BITS;
}

/* the target of the Chebyshev iteration */
typedef struct
{
    ulong R, Q;     /* x* = cos(pi R/Q), 0 < 2R < Q, gcd(R, Q) = 1 */
    int eps;        /* T_Q(x*) = (-1)^R */
    slong kb;       /* bit length of Q */
    slong sb;       /* sin(pi R/Q) > 2^-sb, by Jordan's 2R/Q */
    int r;          /* order */
}
_scp_cheb_struct;

/* limits: the series coefficient words (2j+1) Q +- 1 must fit, j well
   below 2^6 */
#define SCP_CHEB_Q_MAX (UWORD(1) << (FLINT_BITS - 8))
#if FLINT_BITS == 64
#define SCP_DEFAULT_ORDER 12
#define SCP_MAX_ORDER 16
#else
#define SCP_DEFAULT_ORDER 8
#define SCP_MAX_ORDER 10
#endif

/* (v, n + 1) = pi t/q for 0 <= 4t <= q, or its complement pi/2 - pi t/q
   for 4t > q (returning 1): v = pi/4 num/q <= pi/4 with num = 4t resp.
   2(q - 2t), as floor(pi/4 B^n) num / q, below the exact value by less
   than num/q + 1 <= 2 ulps (the product's top limb is below q) */
static int
_scp_v(nn_ptr v, ulong t, ulong q, slong n)
{
    nn_srcptr P;
    ulong num;
    int swap;

    swap = (t > q / 4);     /* 4t > q */
    num = swap ? 2 * (q - 2 * t) : 4 * t;

    P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, n);
    v[n] = mpn_mul_1(v, P, n, num);
    mpn_divrem_1(v, 0, v, n + 1, q);
    FLINT_ASSERT(v[n] == 0);
    return swap;
}

/* the fixed-point kernel path: (ys, n + 1), (yc, n + 1) = sin, cos of
   pi t/q for 0 <= 2t <= q (either output may be NULL) */
static void
_scp_kernel(nn_ptr ys, nn_ptr yc, ulong * err, ulong t, ulong q, slong n)
{
    nn_ptr v, ts, tc;
    ulong e;
    slong n0 = n;
    int swap;
    TMP_INIT;

    /* the kernel takes two limbs at least with 32-bit limbs (as in
       mp_real_sin_cos_bits); the extra limb is truncated at the end */
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif

    TMP_START;
    v = TMP_ALLOC((3 * n + 3) * sizeof(ulong));
    ts = v + n + 1;
    tc = ts + n + 1;

    swap = _scp_v(v, t, q, n);
    _mp_real_sin_cos_kernel(ts, tc, &e, v, n);
    e += 2;

    if (n > n0)
    {
        /* e < B ulps at n limbs are below one at n0, plus the truncation */
        ts += n - n0;
        tc += n - n0;
        e = 2;
    }

    if (swap)
    {
        nn_ptr u = ts;
        ts = tc;
        tc = u;
    }
    if (ys != NULL)
        flint_mpn_copyi(ys, ts, n0 + 1);
    if (yc != NULL)
        flint_mpn_copyi(yc, tc, n0 + 1);
    if (err != NULL)
        *err = e;

    TMP_END;
}

/* res = the ball (y, n + 1) B^-n +- err ulps (err >= 1) */
static void
_scp_fixed_to_ball(mp_real_t res, nn_srcptr y, slong n, ulong err)
{
    mp_real_fit_length(res, n + 1);
    flint_mpn_copyi(res->d, y, n + 1);
    _mp_real_elem_finish(res, n, FLINT_MAX(err, 1), 0);
}

/* the kernel as a ball for cos(pi R/Q), absolute accuracy about 2^-P */
static void
_scp_cheb_kernel(mp_real_t xs, const _scp_cheb_struct * c, slong P)
{
    slong n = (P + 16 + FLINT_BITS - 1) / FLINT_BITS;
    ulong err;
    nn_ptr yc;
    TMP_INIT;

    TMP_START;
    yc = TMP_ALLOC((n + 1) * sizeof(ulong));
    _scp_kernel(NULL, yc, &err, c->R, c->Q, n);
    _scp_fixed_to_ball(xs, yc, n, err);
    TMP_END;
}

/* the most series terms of a step: the common denominators m! (2m-1)!!
   and m! (2m+1)!! must fit a word (about r/2 terms are used) */
#if FLINT_BITS == 64
#define SCP_MAX_TERMS 11
#else
#define SCP_MAX_TERMS 6
#endif

/* D = prod_{j <= i < m} f(i) as a word (0 if it overflows) with
   f(i) = (i+1)(2i+1) (G) resp. (2i+3)(i+1) (F) */
static ulong
_scp_den(slong j, slong m, int F)
{
    ulong D = 1, hi, f;
    slong i;

    for (i = j; i < m; i++)
    {
        f = F ? (ulong) (2 * i + 3) * (ulong) (i + 1) : (ulong) (i + 1) * (ulong) (2 * i + 1);
        umul_ppmm(hi, D, D, f);
        if (hi != 0)
            return 0;
    }
    return D;
}

/* one step from x0 (exact, |acos(x0) - pi R/Q| < pi/(2Q) checked by
   the caller): xs = cos(pi R/Q) to an absolute accuracy of about
   2^-P.  Returns 0 (xs untouched) if the residual is too large for
   the tail bounds or needs too many terms, which the recursion rules
   out. */
static int
_scp_cheb_step(mp_real_t xs, const mp_real_t x0, const _scp_cheb_struct * c, slong P)
{
    mp_real_t V0, V1, v1, t, u, two, A, z, SG, SF, L;
    mp_real_struct zp[SCP_MAX_TERMS + 1];
    slong p, pl, kb = c->kb, i, j, ua, uy, tailG, tailF, cp, mG, mF, m, stop;
    ulong Q = c->Q, DG, DF;
    int ok = 1, odd = (Q & 1);

    p = (P + 4 + FLINT_BITS - 1) / FLINT_BITS + 1;
    /* the ladder's radius grows by about 2 bits per doubling, and x*
       takes A / Q */
    pl = p + (kb + 4 + FLINT_BITS - 1) / FLINT_BITS;

    mp_real_init(V0);
    mp_real_init(V1);
    mp_real_init(v1);
    mp_real_init(t);
    mp_real_init(u);
    mp_real_init(two);
    mp_real_init(A);
    mp_real_init(z);
    mp_real_init(SG);
    mp_real_init(SF);
    mp_real_init(L);
    for (j = 0; j <= SCP_MAX_TERMS; j++)
        mp_real_init(zp + j);

    mp_real_set_ui(two, 2);

    /* (V0, V1) = (V_k, V_(k+1)) for the leading bits k of Q, from k = 1
       (the leading doublings of the short x0 are exact short products);
       for odd Q the ladder stops at k = (Q - 1)/2 */
    stop = odd ? 1 : 0;
    mp_real_mul_2exp_si(v1, x0, 1);
    mp_real_set(V0, v1);
    mp_real_mul(V1, v1, v1, pl);
    mp_real_sub(V1, V1, two, pl);
    for (i = kb - 2; i >= stop; i--)
    {
        mp_real_mul(t, V0, V1, pl);
        mp_real_sub(t, t, v1, pl);
        if ((Q >> i) & 1)
        {
            mp_real_swap(V0, t);
            mp_real_mul(V1, V1, V1, pl);
            mp_real_sub(V1, V1, two, pl);
        }
        else
        {
            mp_real_swap(V1, t);
            mp_real_mul(V0, V0, V0, pl);
            mp_real_sub(V0, V0, two, pl);
        }
    }

    if (odd)
    {
        /* Q = 2k + 1: A = T_k^2 - T_(k+1)^2 = (V_k - V_(k+1))(V_k + V_(k+1))/4
           and V_Q = V_k V_(k+1) - V_1, two products in place of the
           last doubling and x0 T_Q */
        mp_real_sub(u, V0, V1, pl);
        mp_real_add(t, V0, V1, pl);
        mp_real_mul(A, u, t, pl);
        mp_real_mul_2exp_si(A, A, -2);
        mp_real_mul(t, V0, V1, pl);
        mp_real_sub(V0, t, v1, pl);
    }
    else
    {
        /* A = x0 T_Q - T_(Q+1) = (x0 V_Q - V_(Q+1))/2 */
        mp_real_mul(A, x0, V0, pl);
        mp_real_sub(A, A, V1, pl);
        mp_real_mul_2exp_si(A, A, -1);
    }
    /* a = eps A */
    if (c->eps < 0)
        mp_real_neg(A, A);

    /* u = 2 - eps V_Q = 2 (1 - eps C); y = 1 - C^2 = u (4 - u) / 4 and
       z = y / (4 Q^2) = u (4 - u) / (16 Q^2) */
    if (c->eps > 0)
        mp_real_sub(u, two, V0, pl);
    else
        mp_real_add(u, two, V0, pl);
    uy = mp_real_abs_bound_lt_2exp_si(u);
    cp = FLINT_MAX(2, (P + 8 + FLINT_MAX(uy, -P)) / FLINT_BITS + 2);
    mp_real_set_ui(t, 4);
    mp_real_sub(t, t, u, cp);
    mp_real_mul(z, u, t, cp);
    mp_real_mul_2exp_si(z, z, -2);

    uy = mp_real_abs_bound_lt_2exp_si(z);       /* |y| < 2^uy */
    ua = mp_real_abs_bound_lt_2exp_si(A);

    if (uy > -1)
    {
        ok = 0;
        goto cleanup;
    }

    /* the terms: |c_j y^j| < 2^(j uy - 2 kb + 1) for j = 1..mG,
       |(a/Q) f_j y^j| < 2^(ua + j uy - kb + 1) for j = 1..mF, until
       below 2^(-P-4); the tails |y|^(mG+1) / Q^2 resp.
       2 |a| |y|^(mF+1) / Q */
    for (mG = 0; (mG + 1) * uy - 2 * kb + 3 >= -P - 4; mG++)
        ;
    for (mF = 0; ua + (mF + 1) * uy - kb + 1 >= -P - 4; mF++)
        ;
    tailG = (mG + 1) * uy - 2 * kb + 2;
    tailF = ua + (mF + 1) * uy - kb + 2;
    m = FLINT_MAX(mG, mF);

    DG = _scp_den(1, mG, 0);
    DF = _scp_den(0, mF, 1);
    if (m > SCP_MAX_TERMS || DG == 0 || DF == 0)
    {
        ok = 0;
        goto cleanup;
    }

    if (m >= 1)
    {
        /* z = y / (4 Q^2), and its powers at the precision of the terms
           they enter */
        cp = FLINT_MAX(2, (P + 8 + FLINT_MAX(uy - 2 * kb + 3, ua + uy - kb + 1)) / FLINT_BITS + 2);
        mp_real_mul_2exp_si(z, z, -2);
        if (Q <= (UWORD(1) << (FLINT_BITS / 2)) - 1)
            mp_real_div_ui(zp + 1, z, Q * Q, cp);
        else
        {
            mp_real_div_ui(zp + 1, z, Q, cp);
            mp_real_div_ui(zp + 1, zp + 1, Q, cp);
        }
        for (j = 2; j <= m; j++)
        {
            slong e = WORD_MIN;
            if (j <= mG)
                e = j * uy - 2 * kb + 3;
            if (j <= mF)
                e = FLINT_MAX(e, ua + j * uy - kb + 1);
            cp = FLINT_MAX(2, (P + 8 + e) / FLINT_BITS + 2);
            mp_real_mul(zp + j, zp + j - 1, zp + 1, cp);
        }
    }

    /* W1 = G(y) = SG / DG, SG = sum_{j=1}^{mG} Gam_j z^j with the integers
       Gam_j = gamma_j DG, gamma_1 = -2, gamma_(j+1) = gamma_j 2 (4 j^2 Q^2 - 1)
       / ((j+1)(2j+1)): Gam_j = -2 prod_{i<j} 2 (2iQ - 1)(2iQ + 1)
       prod_{j<=i<mG} (i+1)(2i+1) (each term then costs a product by a
       short exact integer) */
    mp_real_zero(SG);
    mp_real_set_ui(L, 2);
    for (j = 1; j <= mG; j++)
    {
        cp = FLINT_MAX(2, (P + 8 + j * uy - 2 * kb + 3) / FLINT_BITS + 2);
        mp_real_mul_ui(t, L, _scp_den(j, mG, 0), SCP_MAX_TERMS * 4 + 4);
        mp_real_mul(t, zp + j, t, cp);
        mp_real_sub(SG, SG, t, FLINT_MAX(2, (P + 8 + uy - 2 * kb + 3) / FLINT_BITS + 2));
        mp_real_mul_ui(L, L, 2 * (2 * j * Q - 1), SCP_MAX_TERMS * 4 + 4);
        mp_real_mul_ui(L, L, 2 * j * Q + 1, SCP_MAX_TERMS * 4 + 4);
    }

    /* W2 = F(y) - 1 = SF / DF, SF = sum_{j=1}^{mF} Phi_j z^j with
       Phi_j = prod_{i<j} 2 ((2i+1) Q - 1)((2i+1) Q + 1)
       prod_{j<=i<mF} (2i+3)(i+1) */
    mp_real_zero(SF);
    mp_real_set_ui(L, 2 * (Q - 1));
    mp_real_mul_ui(L, L, Q + 1, SCP_MAX_TERMS * 4 + 4);
    for (j = 1; j <= mF; j++)
    {
        cp = FLINT_MAX(2, (P + 8 + ua + j * uy - kb + 1) / FLINT_BITS + 2);
        mp_real_mul_ui(t, L, _scp_den(j, mF, 1), SCP_MAX_TERMS * 4 + 4);
        mp_real_mul(t, zp + j, t, cp);
        mp_real_add(SF, SF, t, FLINT_MAX(2, (P + 8 + ua + uy - kb + 1) / FLINT_BITS + 2));
        mp_real_mul_ui(L, L, 2 * ((2 * j + 1) * Q - 1), SCP_MAX_TERMS * 4 + 4);
        mp_real_mul_ui(L, L, (2 * j + 1) * Q + 1, SCP_MAX_TERMS * 4 + 4);
    }

    /* x* = x0 + x0 SG / DG + (a + a SF / DF) / Q, plus the tails */
    if (mF >= 1)
    {
        cp = FLINT_MAX(2, (P + 8 + ua + uy - kb + 1) / FLINT_BITS + 2);
        mp_real_mul(t, A, SF, cp);
        mp_real_div_ui(t, t, DF, cp);
        cp = FLINT_MAX(2, (P + 8 + ua) / FLINT_BITS + 2);
        mp_real_add(t, A, t, cp);
    }
    else
    {
        cp = FLINT_MAX(2, (P + 8 + ua) / FLINT_BITS + 2);
        mp_real_set(t, A);
    }
    mp_real_div_ui(t, t, Q, cp);
    if (mG >= 1)
    {
        cp = FLINT_MAX(2, (P + 8 + uy - 2 * kb + 3) / FLINT_BITS + 2);
        mp_real_mul(u, x0, SG, cp);
        mp_real_div_ui(u, u, DG, cp);
        mp_real_add(t, t, u, p);
    }
    mp_real_add(xs, x0, t, p);
    mp_real_add_error_2exp_si(xs, tailG);
    mp_real_add_error_2exp_si(xs, tailF);

cleanup:
    mp_real_clear(V0);
    mp_real_clear(V1);
    mp_real_clear(v1);
    mp_real_clear(t);
    mp_real_clear(u);
    mp_real_clear(two);
    mp_real_clear(A);
    mp_real_clear(z);
    mp_real_clear(SG);
    mp_real_clear(SF);
    mp_real_clear(L);
    for (j = 0; j <= SCP_MAX_TERMS; j++)
        mp_real_clear(zp + j);

    return ok;
}

/* whether a Chebyshev step (and the levels below it) beats the kernel
   at P bits */
static int
_scp_cheb_worth(slong P, slong kb)
{
    return kb <= _scp_cheb_kbmax((P + FLINT_BITS - 1) / FLINT_BITS);
}

/* xs = cos(pi R/Q) to an absolute accuracy of about 2^-P (a ball),
   by the kernel or by a step from the recursion */
static void
_scp_cheb_rec(mp_real_t xs, const _scp_cheb_struct * c, slong P, int force)
{
    slong need, k, erho;
    int r = c->r, ok = 0;
    mp_real_t xb, x0;

    /* the accuracy of the approximation: the step leaves about
       r (need - kb - sb - 2) bits; enough for the branch check and a
       small residual */
    need = (P + 16 + (r - 1) * (c->kb + c->sb + 2) + r - 1) / r;
    need = FLINT_MAX(need, c->kb + 2 * c->sb + 32);

    if (need >= P - 32 || !(force || _scp_cheb_worth(P, c->kb)))
    {
        _scp_cheb_kernel(xs, c, P);
        return;
    }

    mp_real_init(xb);
    mp_real_init(x0);

    _scp_cheb_rec(xb, c, need, 0);

    /* x0 = the midpoint of xb truncated to k fraction limbs (xb
       encloses x* in [0, 1), so its exponent is at most 1) */
    k = (need + 8 + FLINT_BITS - 1) / FLINT_BITS + 1;
    mp_real_set(x0, xb);
    x0->err = 0;
    erho = -FLINT_BITS * k;
    if (x0->size > 0 && x0->exp + k < x0->size)
    {
        slong keep = x0->exp + k, lo;

        if (keep <= 0)
            mp_real_zero(x0);
        else
        {
            flint_mpn_copyi(x0->d, x0->d + x0->size - keep, keep);
            x0->size = keep;
            for (lo = 0; x0->d[lo] == 0; lo++)
                ;
            if (lo > 0)
            {
                flint_mpn_copyi(x0->d, x0->d + lo, keep - lo);
                x0->size = keep - lo;
            }
        }
    }

    /* |x0 - x*| < rho = rad(xb) + B^-k < 2^erho; the branch condition
       |acos(x0) - pi R/Q| < pi/(2Q) holds when rho <= min(s^2/8, s/Q)
       with s = sin(pi R/Q) > 2^-sb (then |x| + rho <= 1 - 3 s^2/8 and
       sqrt(1 - x^2) >= 0.78 s on the interval) */
    if (xb->err != 0)
        erho = FLINT_MAX(erho, FLINT_BITS * (xb->exp - xb->size)
            + (slong) FLINT_BIT_COUNT(xb->err));
    erho += 1;

    if (x0->size != 0 && !x0->negative && x0->exp <= 0
        && erho <= -(2 * c->sb + 3) && erho <= -(c->sb + c->kb))
        ok = _scp_cheb_step(xs, x0, c, P);

    if (!ok)
        _scp_cheb_kernel(xs, c, P);

    mp_real_clear(xb);
    mp_real_clear(x0);
}

static void
_scp_cheb_init(_scp_cheb_struct * c, ulong R, ulong Q, int r)
{
    c->R = R;
    c->Q = Q;
    c->eps = (R & 1) ? -1 : 1;
    c->kb = FLINT_BIT_COUNT(Q);
    c->sb = FLINT_BIT_COUNT(Q / (2 * R));
    c->r = (r <= 0) ? SCP_DEFAULT_ORDER : FLINT_MAX(FLINT_MIN(r, SCP_MAX_ORDER), 2);
}

/* balls cs = cos(pi t/q), sn = sin(pi t/q) (either may be NULL) to an
   absolute accuracy of about B^-n, 0 < 2t < q, gcd(t, q) = 1, q > 6,
   by the Chebyshev iteration on cos(pi t/q) (target A) or on
   sin(pi t/q) = cos(pi (q - 2t)/(2q)) (target B), whichever costs
   less: B has the larger denominator for odd q but is well conditioned
   for t/q < 1/4 */
static void
_scp_cheb_ball(mp_real_t cs, mp_real_t sn, ulong t, ulong q, slong n, int r)
{
    _scp_cheb_struct cA, cB, * c;
    mp_real_t xs, ss;
    slong P;
    int useB, need_s;

    _scp_cheb_init(&cA, t, q, r);

    useB = 0;
    if (q <= SCP_CHEB_Q_MAX / 2)
    {
        ulong R = q - 2 * t, Q = 2 * q, g = n_gcd(R, Q);
        double costA, costB;

        _scp_cheb_init(&cB, R / g, Q / g, r);
        /* the ladder's bits times the working precision (the guard: the
           ladder, and the conditioning of the square root) */
        costA = (double) cA.kb * (FLINT_BITS * n + cA.kb + cA.sb + 16);
        costB = (double) cB.kb * (FLINT_BITS * n + cB.kb + cB.sb + 16);
        useB = (costB <= costA);
    }
    c = useB ? &cB : &cA;

    /* x = cos phi* is cos for A and sin for B, s = sin phi* the other
       way */
    need_s = useB ? (cs != NULL) : (sn != NULL);

    mp_real_init(xs);
    mp_real_init(ss);

    /* the square root amplifies the error of x by x/s < 2^sb (and
       1 - x^2 >= s^2 must be accurate enough for mp_real_sqrt) */
    P = FLINT_BITS * n + 8 + (need_s ? 2 * c->sb + 32 : 0);
    _scp_cheb_rec(xs, c, P, 1);

    if (need_s)
    {
        mp_real_t one;
        slong p = n + 2 + (c->sb + FLINT_BITS - 1) / FLINT_BITS;

        mp_real_init(one);
        mp_real_set_ui(one, 1);
        mp_real_mul(ss, xs, xs, p);
        mp_real_sub(ss, one, ss, p);
        mp_real_sqrt(ss, ss, p);
        mp_real_clear(one);
    }

    if (useB)
    {
        if (sn != NULL)
            mp_real_swap(sn, xs);
        if (cs != NULL)
            mp_real_swap(cs, ss);
    }
    else
    {
        if (cs != NULL)
            mp_real_swap(cs, xs);
        if (sn != NULL)
            mp_real_swap(sn, ss);
    }

    mp_real_clear(xs);
    mp_real_clear(ss);
}

/* the same as fixed-point numbers (ys, n + 1), (yc, n + 1) */
static void
_scp_cheb(nn_ptr ys, nn_ptr yc, ulong * err, ulong t, ulong q, slong n, int r)
{
    mp_real_t cs, sn;
    ulong e1 = 0, e2 = 0;

    mp_real_init(cs);
    mp_real_init(sn);
    _scp_cheb_ball(yc ? cs : NULL, ys ? sn : NULL, t, q, n, r);
    if (ys != NULL)
        e1 = _mp_real_ball_get_fixed1(ys, sn, n);
    if (yc != NULL)
        e2 = _mp_real_ball_get_fixed1(yc, cs, n);
    if (err != NULL)
        *err = FLINT_MAX(e1, e2);
    mp_real_clear(cs);
    mp_real_clear(sn);
}

/* (res, n + 1) = the constant c (units limb included), exact */
static void
_scp_set_fixed_ui(nn_ptr res, slong n, ulong c)
{
    if (res != NULL)
    {
        flint_mpn_zero(res, n);
        res[n] = c;
    }
}

/* (res, n + 1) = 1/2 */
static void
_scp_set_fixed_half(nn_ptr res, slong n)
{
    if (res != NULL)
    {
        flint_mpn_zero(res, n + 1);
        res[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    }
}

/* (res, n + 1) = sqrt(2)/2 = rsqrt(2) resp. sqrt(3)/2 = 3 rsqrt(3) / 2
   within 4 ulps (the root within 2, times 3/2, and the shift) */
static void
_scp_set_fixed_sqrt(nn_ptr res, slong n, ulong a)
{
    if (res != NULL)
    {
        _mp_real_rsqrt_ui_newton(res, a, n);
        res[n] = 0;
        if (a == 3)
        {
            res[n] = mpn_mul_1(res, res, n, 3);
            mpn_rshift(res, res, n + 1, 1);
        }
    }
}

/* the closed forms below beat the kernel from about 12 limbs (24 for
   q = 8, with two square roots; measured on 64-bit x86, in bits here),
   and the Chebyshev iteration by 1.4 to 5 times */
#define SCP_QUADRATIC_MIN_BITS(q) (((q) == 8) ? 24 * 64 : 12 * 64)

/* sin and cos of pi t/q for q = 5, 8, 10, 12 (0 < 2t < q, gcd(t, q) = 1)
   by square roots, in ball arithmetic: with r5 = sqrt(5), r2 = sqrt(2),
   r6 = sqrt(6),

       cos(pi/5) = sin(3 pi/10) = (r5 + 1)/4,
       cos(2 pi/5) = sin(pi/10) = (r5 - 1)/4,
       sin(pi/5) = cos(3 pi/10) = sqrt((5 - r5)/8),
       sin(2 pi/5) = cos(pi/10) = sqrt((5 + r5)/8),
       sin(pi/12), cos(pi/12) = (r6 -+ r2)/4,
       cos(pi/8) = sqrt(2 + r2)/2 = (2 + r2) rho/2,
       sin(pi/8) = sqrt(2 - r2)/2 = r2 rho/2,   rho = 1/sqrt(2 + r2),

   the other t swapping sin and cos */
static void
_scp_quadratic(nn_ptr ys, nn_ptr yc, ulong * err, ulong t, ulong q, slong n)
{
    mp_real_t lo, hi, u, v;
    slong w = n + 2;
    ulong e1 = 0, e2 = 0;
    int swap;
    nn_ptr ylo, yhi;

    mp_real_init(lo);
    mp_real_init(hi);
    mp_real_init(u);
    mp_real_init(v);

    /* lo = the smaller of the values, hi = the larger; for q = 5, 10 the
       rational one and the square root, by the target */
    if (q == 5 || q == 10)
    {
        /* rat = (r5 + e)/4 and root = sqrt((5 - e r5)/8), e = +1 when
           the rational value is the larger (cos(pi/5), sin(3 pi/10)) */
        int e = (q == 5) ? (t == 1) : (t == 3);
        int want_rat = (q == 5) ? (yc != NULL) : (ys != NULL);
        int want_root = (q == 5) ? (ys != NULL) : (yc != NULL);

        mp_real_rsqrt_ui(u, 5, w);
        mp_real_mul_ui(u, u, 5, w);
        if (want_rat)
        {
            mp_real_set_ui(v, 1);
            if (e)
                mp_real_add(lo, u, v, w);
            else
                mp_real_sub(lo, u, v, w);
            mp_real_mul_2exp_si(lo, lo, -2);
        }
        if (want_root)
        {
            mp_real_set_ui(v, 5);
            if (e)
                mp_real_sub(hi, v, u, w);
            else
                mp_real_add(hi, v, u, w);
            mp_real_mul_2exp_si(hi, hi, -3);
            mp_real_sqrt(hi, hi, w);
        }
        /* lo holds the rational value, hi the root: cos, sin for q = 5 */
        swap = (q == 10);
    }
    else if (q == 12)
    {
        mp_real_rsqrt_ui(u, 6, w);
        mp_real_mul_ui(u, u, 6, w);
        mp_real_rsqrt_ui(v, 2, w);
        mp_real_mul_2exp_si(v, v, 1);
        mp_real_sub(lo, u, v, w);
        mp_real_add(hi, u, v, w);
        mp_real_mul_2exp_si(lo, lo, -2);
        mp_real_mul_2exp_si(hi, hi, -2);
        /* hi = cos, lo = sin for t = 1 */
        swap = (t == 1);
    }
    else
    {
        FLINT_ASSERT(q == 8);
        mp_real_rsqrt_ui(v, 2, w);
        mp_real_mul_2exp_si(v, v, 1);
        mp_real_set_ui(u, 2);
        mp_real_add(u, u, v, w);
        mp_real_rsqrt(hi, u, w);
        mp_real_mul(lo, v, hi, w);
        mp_real_mul(hi, u, hi, w);
        mp_real_mul_2exp_si(lo, lo, -1);
        mp_real_mul_2exp_si(hi, hi, -1);
        /* hi = cos, lo = sin for t = 1 */
        swap = (t == 1);
    }

    /* unswapped: lo -> yc, hi -> ys */
    ylo = swap ? ys : yc;
    yhi = swap ? yc : ys;
    if (ylo != NULL)
        e1 = _mp_real_ball_get_fixed1(ylo, lo, n);
    if (yhi != NULL)
        e2 = _mp_real_ball_get_fixed1(yhi, hi, n);
    if (err != NULL)
        *err = FLINT_MAX(e1, e2);

    mp_real_clear(lo);
    mp_real_clear(hi);
    mp_real_clear(u);
    mp_real_clear(v);
}

static void _scp_lowest(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q,
    slong n, int alg, int r);

void
_mp_real_sin_cos_pi_ui_div_ui_tune(nn_ptr ys, nn_ptr yc, ulong * err,
    ulong p, ulong q, slong n, int alg, int r)
{
    ulong g;

    FLINT_ASSERT(q >= 1);
    FLINT_ASSERT(p <= q / 2);
    FLINT_ASSERT(n >= 1);

    if (ys == NULL && yc == NULL)
        return;

    if (p == 0)
    {
        _scp_set_fixed_ui(ys, n, 0);
        _scp_set_fixed_ui(yc, n, 1);
        if (err != NULL)
            *err = 0;
        return;
    }

    g = n_gcd(p, q);
    _scp_lowest(ys, yc, err, p / g, q / g, n, alg, r);
}

/* the dispatch for a fraction in lowest terms with 0 < p <= q/2 */
static void
_scp_lowest(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n,
    int alg, int r)
{
    ulong e;

    /* exact and quadratic values */
    if (q == 2)
    {
        _scp_set_fixed_ui(ys, n, 1);
        _scp_set_fixed_ui(yc, n, 0);
        e = 0;
    }
    else if (q == 3)
    {
        _scp_set_fixed_sqrt(ys, n, 3);
        _scp_set_fixed_half(yc, n);
        e = 4;
    }
    else if (q == 4)
    {
        _scp_set_fixed_sqrt(ys, n, 2);
        _scp_set_fixed_sqrt(yc, n, 2);
        e = 4;
    }
    else if (q == 6)
    {
        _scp_set_fixed_half(ys, n);
        _scp_set_fixed_sqrt(yc, n, 3);
        e = 4;
    }
    else if ((q == 5 || q == 8 || q == 10 || q == 12) && (alg == 3
            || (alg == 0 && FLINT_BITS * n >= SCP_QUADRATIC_MIN_BITS(q))))
    {
        _scp_quadratic(ys, yc, &e, p, q, n);
    }
    else
    {
        if (alg == 0)
            alg = (q <= SCP_CHEB_Q_MAX
                && (slong) FLINT_BIT_COUNT(q) <= _scp_cheb_kbmax(n)) ? 2 : 1;

        if (alg == 2 && q <= SCP_CHEB_Q_MAX)
            _scp_cheb(ys, yc, &e, p, q, n, r);
        else
            _scp_kernel(ys, yc, &e, p, q, n);
    }

    if (err != NULL)
        *err = e;
}

void
_mp_real_sin_cos_pi_ui_div_ui(nn_ptr ys, nn_ptr yc, ulong * err,
    ulong p, ulong q, slong n)
{
    _mp_real_sin_cos_pi_ui_div_ui_tune(ys, yc, err, p, q, n, 0, 0);
}

/* res = (-1)^neg c for c = 0, 1 or 1/2 */
static void
_scp_set_exact(mp_real_t res, int c2, int neg)
{
    /* c2 = 2c */
    if (c2 == 0)
        mp_real_zero(res);
    else
    {
        mp_real_set_ui(res, 1);
        if (c2 == 1)
            mp_real_mul_2exp_si(res, res, -1);
        if (neg)
            mp_real_neg(res, res);
    }
}

void
mp_real_sin_cos_pi_ui_div_ui(mp_real_t res1, mp_real_t res2, ulong p,
    ulong q, slong n)
{
    ulong g, a, t, err;
    int negs, negc, exact1, exact2;
    slong w;

    if (q == 0)
        flint_throw(FLINT_ERROR, "mp_real_sin_cos_pi_ui_div_ui: q = 0\n");

    n = FLINT_MAX(n, 1);

    /* lowest terms (0/1 for p = 0) */
    g = n_gcd(p, q);
    p /= g;
    q /= g;

    /* pi p/q = pi a + pi t/q, 0 <= t < q; then pi t/q = pi - pi t'/q */
    a = p / q;
    t = p % q;
    negs = negc = (int) (a & 1);
    if (t > q - t)
    {
        t = q - t;
        negc ^= 1;
    }

    /* the rational values 0, 1 and 1/2 directly: sin for q = 1, 2, 6,
       cos for q = 1, 2, 3 */
    exact1 = (q <= 2 || q == 6);
    exact2 = (q <= 2 || q == 3);
    if (res1 != NULL && exact1)
    {
        _scp_set_exact(res1, (q == 1) ? 0 : (q == 2) ? 2 : 1, negs);
        res1 = NULL;
    }
    if (res2 != NULL && exact2)
    {
        _scp_set_exact(res2, (q == 2) ? 0 : (q == 1) ? 2 : 1, negc);
        res2 = NULL;
    }
    if (res1 == NULL && res2 == NULL)
        return;

    /* the working precision: |sin| >= 2t/q, |cos| >= (q - 2t)/q, both
       at least 1/q here (nonzero), and the fixed-point bound is at most
       13 bits, so bits(q) + 16 extra bits give the relative accuracy
       B^-n */
    w = n + (FLINT_BIT_COUNT(q) + 16 + FLINT_BITS - 1) / FLINT_BITS;

    if (res1 != NULL)
        mp_real_fit_length(res1, w + 1);
    if (res2 != NULL)
        mp_real_fit_length(res2, w + 1);

    _scp_lowest(res1 ? res1->d : NULL, res2 ? res2->d : NULL, &err, t, q, w, 0, 0);

    if (res1 != NULL)
        _mp_real_elem_finish(res1, w, FLINT_MAX(err, 1), negs);
    if (res2 != NULL)
        _mp_real_elem_finish(res2, w, FLINT_MAX(err, 1), negc);
}

/* tangent ********************************************************************/

/* T = tan(pi t/q) for 0 < 2t < q, gcd(t, q) = 1, q >= 3, q != 4, to a
   relative accuracy of about B^-w (w >= 2 + the bits of q in limbs):
   closed forms for q = 3, 6, 8, 12 and (from the size at which they
   beat the kernel for sin and cos) q = 5, 10; the Chebyshev sine and
   cosine divided, where they beat the kernel; else the tan kernel on
   v = pi t/q, or on its complement pi/2 - v for 4t > q, then inverted */
static void
_scp_tan_ball(mp_real_t T, ulong t, ulong q, slong w)
{
    if (q == 3 || q == 6 || q == 8 || q == 12)
    {
        /* tan(pi/3) = sqrt(3), tan(pi/6) = 1/sqrt(3), tan(pi/8), tan(3 pi/8)
           = sqrt(2) -+ 1, tan(pi/12), tan(5 pi/12) = 2 -+ sqrt(3) */
        mp_real_t u;
        mp_real_init(u);
        mp_real_rsqrt_ui(T, (q == 8) ? 2 : 3, w);
        if (q == 3)
            mp_real_mul_ui(T, T, 3, w);
        else if (q == 8)
        {
            mp_real_mul_2exp_si(T, T, 1);
            mp_real_set_ui(u, 1);
            if (t == 1)
                mp_real_sub(T, T, u, w);
            else
                mp_real_add(T, T, u, w);
        }
        else if (q == 12)
        {
            mp_real_mul_ui(T, T, 3, w);
            mp_real_set_ui(u, 2);
            if (t == 1)
                mp_real_sub(T, u, T, w);
            else
                mp_real_add(T, u, T, w);
        }
        mp_real_clear(u);
    }
    else if ((q == 5 || q == 10) && FLINT_BITS * w >= SCP_QUADRATIC_MIN_BITS(q))
    {
        /* tan(pi/5), tan(2 pi/5) = sqrt(5 -+ 2 sqrt(5)),
           tan(pi/10), tan(3 pi/10) = sqrt(1 -+ 2/sqrt(5)) */
        mp_real_t u;
        int minus = (t == 1);

        mp_real_init(u);
        mp_real_rsqrt_ui(T, 5, w);
        mp_real_mul_2exp_si(T, T, 1);
        if (q == 5)
        {
            mp_real_mul_ui(T, T, 5, w);
            mp_real_set_ui(u, 5);
        }
        else
            mp_real_set_ui(u, 1);
        if (minus)
            mp_real_sub(T, u, T, w);
        else
            mp_real_add(T, u, T, w);
        mp_real_sqrt(T, T, w);
        mp_real_clear(u);
    }
    else if (q <= SCP_CHEB_Q_MAX && (slong) FLINT_BIT_COUNT(q) <= _scp_cheb_kbmax(w))
    {
        mp_real_t cs;
        mp_real_init(cs);
        _scp_cheb_ball(cs, T, t, q, w, 0);
        mp_real_div(T, T, cs, w);
        mp_real_clear(cs);
    }
    else
    {
        nn_ptr v, y;
        ulong err;
        int swap;
        TMP_INIT;

        TMP_START;
        v = TMP_ALLOC((w + 1) * sizeof(ulong));
        y = TMP_ALLOC((w + 1) * sizeof(ulong));
        swap = _scp_v(v, t, q, w);

        /* the kernel within its bound, v's 2 ulps times sec^2 v <= 2 */
        _mp_real_tan_kernel(y, &err, v, w);
        if (!swap)
        {
            mp_real_fit_length(T, w + 1);
            flint_mpn_copyi(T->d, y, w + 1);
            _mp_real_elem_finish(T, w, err + 5, 0);
        }
        else if (!_mp_real_inv_mpn(T, y, w + 1, -w, err + 5, WORD_MIN, w, 0))
            flint_throw(FLINT_ERROR, "mp_real_tan_pi_ui_div_ui: tan v not positive\n");
        TMP_END;
    }
}

int
mp_real_tan_pi_ui_div_ui(mp_real_t res, ulong p, ulong q, slong n)
{
    ulong g, t;
    slong w;
    int neg = 0;
    mp_real_t T;

    if (q == 0)
        flint_throw(FLINT_ERROR, "mp_real_tan_pi_ui_div_ui: q = 0\n");

    n = FLINT_MAX(n, 1);

    /* lowest terms; tan has period pi: tan(pi p/q) = tan(pi t/q) */
    g = n_gcd(p, q);
    p /= g;
    q /= g;
    t = p % q;

    if (t == 0)
    {
        mp_real_zero(res);
        return 1;
    }

    /* a pole at pi/2 */
    if (q == 2)
    {
        mp_real_zero(res);
        return 0;
    }

    /* tan(pi t/q) = -tan(pi (q - t)/q), now 0 < 2t < q */
    if (t > q - t)
    {
        t = q - t;
        neg = 1;
    }

    if (q == 4)
    {
        mp_real_set_ui(res, 1);
        if (neg)
            mp_real_neg(res, res);
        return 1;
    }

    /* the value is at least tan(pi/q) > pi/q, and at most its reciprocal,
       so an absolute B^-w (with the bits of q and 16 more in guard) is
       a relative B^-n */
    w = n + 1 + (FLINT_BIT_COUNT(q) + 16 + FLINT_BITS - 1) / FLINT_BITS;

    mp_real_init(T);
    _scp_tan_ball(T, t, q, w);
    if (neg)
        mp_real_neg(T, T);
    mp_real_swap(res, T);
    mp_real_clear(T);
    return 1;
}
