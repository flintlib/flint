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
#include "fixed.h"
#include "impl.h"

#include "tan_rotate.inc"

/* the sizes for which the rotation is emitted in registers */
#define FIXED_TAN_ROTATE_NMAX 7

/* The tangent half-angle reconstruction.

   The angle is reduced exactly as for fixed_sin_cos_bitwise_rs: with the
   cached angles A_i = atan(2^-i),

       x/2 = sum_{i in used} A_i + t',        t' < 2^-r

   and the rotation W = prod (1 + i 2^-i) over the used indices is built
   in the same way, two shifts and two add/subtracts per factor.

   The difference is what is done with it.  Since arg(W) = sum A_i,

       tan(sum A_i) = Im W / Re W = wy / wx,

   and because TAN IS A RATIO, the growth of |W| cancels: no |W|^2, no
   normalizing square root, no unimodular correction -- the whole reason
   for the conjugate-ratio identity used by sin/cos simply evaporates.
   With u = tan(t') the addition formula gives the half-angle tangent in
   ONE division,

       t = tan(x/2) = (wy + wx u) / (wx - wy u),

   whose denominator is close to wx in (0.72, 1.17): no cancellation.
   Everything else follows from t by the half-angle formulas

       sin x = 2t / (1 + t^2),
       cos x = (1 - t^2) / (1 + t^2),
       tan x = 2t / (1 - t^2),

   sharing the single squaring t^2.  Note 0 <= t < tan(1/2) = 0.5464, so
   2t may EXCEED one and the outputs carry a unit limb.  The divisions
   are carried out with numerator and denominator halved, sin = t/D and
   cos = ((1 - t^2)/2)/D with D = (1 + t^2)/2 in [0.5, 0.6492), and
   tan as twice t/(1 - t^2): every divisor is then a NORMALIZED n-limb
   value (top bit set), one limb shorter than 1 + t^2 itself, and every
   division is well conditioned.

   Measured against the conjugate-ratio reconstruction this is a third to
   a half faster across n = 1..8, and it yields sin, cos or tan
   separately rather than computing two functions and discarding one. */

#if FLINT_BITS == 64

/* the reduction parameter each series was built for, re-swept every
   time the tail got cheaper -- and it did so twice more since the
   original 5, 6, 15, 27, 27, 30: after the register rotation and
   normalized divisors, and again after the one-reciprocal tail; the
   optimum keeps moving down.  The candidates are validated at 8-10k
   adversarial points before being locked: the sweep's 300-point error
   estimate MISSES spikes (r = 17 at n = 5 and r = 15 at n = 7 measured
   fastest but blow the tan budget at scale). */

#include "hand_mulhi.inc"

#else

/* no hand-written tangent series off 64-bit limbs: every size then
   takes tan(t') from the sine and cosine series with one division.
   The register rotation is portable and is used regardless. */
#define FIXED_TAN_NMAX 0

#endif

#if FLINT_BITS == 64

#endif

void
_fixed_tan_halfangle_mid(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r, void (*series)(nn_ptr, nn_srcptr))
{
    slong wn = n + 1, nc, num, i, j;
    nn_ptr t, u, wx, wy, va, vb, T, D, N, Q, sc, ss, cc;
    slong * used;
    TMP_INIT;

    _fixed_atans_ensure(n, r);
    nc = _fixed_atans_n;

    TMP_START;
    t = TMP_ALLOC((n + n + wn + wn + n + n + wn + wn
        + (2 * n + 2) + (2 * n + 2) + (2 * n + 4)
        + (n + 2) + (n + 2)) * sizeof(ulong));
    u = t + n;
    wx = u + n;
    wy = wx + wn;
    va = wy + wn;
    vb = va + n;
    T = vb + n;                 /* wn: t^2 / products */
    D = T + wn;                 /* wn: the divisor */
    N = D + wn;                 /* 2n+2: division numerator */
    Q = N + (2 * n + 2);        /* 2n+2: quotient */
    sc = Q + (2 * n + 2);       /* 2n+4: division scratch */
    ss = sc + (2 * n + 4);      /* n+2 each: sin, cos of the residual */
    cc = ss + (n + 2);
    used = TMP_ALLOC(FIXED_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    /* reduce x/2; the residual t' is wanted as it stands, so unlike the
       sin/cos reduction there is no doubling and no extra index */
    mpn_rshift(t, x, n, 1);
    num = _fixed_bitwise_reduce(t, n, r, 1, _fixed_atans, nc, used);

    /* W = prod (1 + i 2^-i).  For the small sizes this runs in
       registers (tan_rotate.inc): the generic loop below materializes
       wx >> i and wy >> i into scratch each time round -- zeroing them,
       filling them with a shift or a copy, then an mpn_sub and an
       mpn_add -- none of which is necessary when the components fit in
       registers and the shift is a pair of funnel shifts. */
    if (n <= FIXED_TAN_ROTATE_NMAX)
    {
        _fixed_tan_rotate_tab[n](wx, wy, used, num);
    }
    else
    {
        flint_mpn_zero(wx, n);
        wx[n] = 1;
        flint_mpn_zero(wy, wn);

        for (j = 0; j < num; j++)
        {
            slong ii = used[j], q = ii / FLINT_BITS;
            int b = (int) (ii - q * FLINT_BITS);

            flint_mpn_zero(va + (wn - q), q);
            flint_mpn_zero(vb + (n - q), q);
            if (b != 0)
            {
                mpn_rshift(va, wx + q, wn - q, b);
                if (n - q > 0)
                    mpn_rshift(vb, wy + q, n - q, b);
            }
            else
            {
                flint_mpn_copyi(va, wx + q, wn - q);
                if (n - q > 0)
                    flint_mpn_copyi(vb, wy + q, n - q);
            }
            mpn_sub(wx, wx, wn, vb, n);
            mpn_add(wy, wy, wn, va, n);
        }
    }

    /* TT and DE, the imaginary and real parts of W (cos t' + i sin t')
       up to positive real factors that CANCEL in every ratio below just
       as |W| does: TT/DE is tan(x/2).  TT, at most
       tan(1/2) = 0.5464 times DE in (0.72, 1.17), stays below 0.64 and
       never uses its unit limb.

       Where a tangent series exists, u = tan(t') comes from it
       directly and TT = wy + wx u, DE = wx - wy u (dropping the
       common factor cos t'); wx and wy are (n+1)-limb values scaled
       by 2^-64n and u an n-limb fraction, so each product at that
       scale is the top n+1 limbs of the (2n+1)-limb integer product.

       Beyond those sizes the tangent u itself is NEVER DIVIDED OUT:
       writing cos t' = 1 - g with g < 2^(-2r-1),

           TT = wy cos t' + wx sin t' = wy + wx ss - wy g,
           DE = wx cos t' - wy sin t' = wx - wx g - wy ss,

       four mulhighs against the tiny ss and g in place of the division
       ss/cc that used to make u (this was the last of the three
       divisions the old reconstruction chain spent before reaching the
       output reciprocals).  No difference can go negative: TT's
       subtrahend wy g is below wy, DE's two are below wx. */
    if (series != NULL)
    {
        series(u, t);                       /* u = tan(t') */
        flint_mpn_mul(N, wx, wn, u, n);
        mpn_add_n(T, wy, N + n, wn);        /* TT (T[n] = 0) */
        flint_mpn_mul(N, wy, wn, u, n);
        mpn_sub_n(D, wx, N + n, wn);        /* DE */
    }
    else
    {
        fixed_sin_cos_rs(ss, cc, t, n);

        if (cc[n])
            flint_mpn_zero(cc, n);          /* cos t' = 1: g = 0 */
        else
            mpn_neg(cc, cc, n);             /* g = 1 - cos t' */

        flint_mpn_mulhigh_n(va, wx, ss, n); /* wx ss */
        if (wx[n])
            mpn_add_n(va, va, ss, n);       /* < 2^(1-r): no carry */
        flint_mpn_mulhigh_n(vb, wy, cc, n); /* wy g */
        T[n] = mpn_add_n(T, wy, va, n);     /* wy + wx ss < 1 */
        mpn_sub_n(T, T, vb, n);             /* TT */

        flint_mpn_mulhigh_n(va, wx, cc, n); /* wx g */
        if (wx[n])
            mpn_add_n(va, va, cc, n);
        flint_mpn_mulhigh_n(vb, wy, ss, n); /* wy ss */
        flint_mpn_copyi(D, wx, wn);
        mpn_sub(D, D, wn, va, n);
        mpn_sub(D, D, wn, vb, n);           /* DE */
    }

    if (D[n])
    {
        mpn_rshift(D, D, wn, 1);
        mpn_rshift(T, T, n, 1);
    }

    /* The half-angle tangent t = TT/DE is never materialized.  With

           A = TT DE,  B = DE^2,  C = TT^2,  S = B + C = DE^2 (1 + t^2)

       the outputs are

           sin = 2A/S,   cos = 1 - 2C/S,   tan = 2A/(B - C),

       each ONE reciprocal (against a divisor normalized to n limbs
       with the top bit set by a shift of at most two bits, whose
       exponent moves into the final doubling shift) and one mulhigh
       per output -- where computing t first would have spent a whole
       division before the reciprocal.  Every former special case
       dissolves: t^2 = 0 makes C = 0 and cos exactly one through the
       generic code path.

       Ranges after the joint halving: A < 0.55, B in [0.25, 1),
       C < 0.30, S in [0.25, 1.30) with S >= 0.5 whenever the halving
       did not fire (B > 0.518 then), W = B - C = DE^2 (1 - t^2) in
       (0.175, 1).  The doubling shifts amplify the last mulhigh's
       truncation by at most 8 (sin/cos) resp. 16 (tan) ulp, absorbed
       many times over by the 6r + 128 and 8r + 256 budgets. */
    {
        nn_ptr A = u, B = va, C = vb, S = T, R = Q, h = sc, W = D;
        ulong topbit = UWORD(1) << (FLINT_BITS - 1);
        int e;

        if (ysin != NULL || ytan != NULL)
            flint_mpn_mulhigh_n(A, T, D, n);
        flint_mpn_sqrhigh(B, D, n);
        flint_mpn_sqrhigh(C, T, n);

        flint_mpn_store(N, 2 * n, ~UWORD(0));
        N[2 * n - 1] = ~UWORD(0) >> 1;          /* 2^(128n-1) - 1 */

        if (ysin != NULL || ycos != NULL)
        {
            /* S, normalized to [0.5, 1) by a shift of e in {-1,0,1} */
            if (mpn_add_n(S, B, C, n))
            {
                mpn_rshift(S, S, n, 1);
                S[n - 1] |= topbit;
                e = 1;
            }
            else if (!(S[n - 1] & topbit))
            {
                mpn_lshift(S, S, n, 1);
                e = -1;
            }
            else
                e = 0;

            /* R = 1/(2 S'): the all-ones numerator keeps R below
               2^(64n) even when S' = 1/2 exactly */
            mpn_tdiv_qr(R, h, 0, N, 2 * n, S, n);

            if (ysin != NULL)
            {
                /* sin = 2A/S = 2^(2-e) A R < 0.842 */
                flint_mpn_mulhigh_n(ysin, A, R, n);
                ysin[n] = 0;
                mpn_lshift(ysin, ysin, n + 1, 2 - e);
            }
            if (ycos != NULL)
            {
                /* cos = 1 - 2C/S with 2C/S = 2t^2/(1 + t^2) < 0.46 */
                flint_mpn_mulhigh_n(h, C, R, n);
                mpn_lshift(h, h, n, 2 - e);
                if (flint_mpn_zero_p(h, n))
                {
                    flint_mpn_zero(ycos, n);
                    ycos[n] = 1;
                }
                else
                {
                    mpn_neg(ycos, h, n);
                    ycos[n] = 0;
                }
            }
        }

        if (ytan != NULL)
        {
            /* W = B - C > 0.175: at most two leading zero bits */
            mpn_sub_n(W, B, C, n);
            e = flint_clz(W[n - 1]);
            if (e)
                mpn_lshift(W, W, n, e);

            mpn_tdiv_qr(R, h, 0, N, 2 * n, W, n);

            /* tan = 2A/(B - C) = 2^(2+e) A R < 1.56 */
            flint_mpn_mulhigh_n(ytan, A, R, n);
            ytan[n] = 0;
            mpn_lshift(ytan, ytan, n + 1, 2 + e);
        }
    }

    (void) i;
    TMP_END;
    return;
}

/* smallest n at which each reduction parameter becomes optimal for
   the shared half-angle path of fixed_sin_cos_bitwise_rs and
   fixed_tan_bitwise_rs (generated by src/fixed/tune/tune-bitwise-r.c
   with nmax deep enough to reach the r = 768 end of the ladder;
   x86-64 defaults) */
static const int _fixed_trig_bitwise_rs_r_tab[] =
    {32, 64, 128, 192, 256, 320, 384, 448, 512, 576, 640, 704, 768};
static const short _fixed_trig_bitwise_rs_n_tab[] =
    {1, 11, 40, 59, 254, 301, 320, 919, 920, 1198, 1437, 1449, 1450};

int
fixed_trig_bitwise_rs_default_r(slong n)
{
    slong j;

#if FLINT_BITS == 64
    /* the compile-time constants of trig_opt_<n>.c: keep in sync
       with the dev/tune_fixed.py --pin values that emitted them */
    static const int opt_r[] =
        {0, 4, 5, 9, 14, 15, 18, 16, 16, 16, 19, 23, 25};
    if (n <= 12)
        return opt_r[n];
#endif

    for (j = 0; j + 1 < (slong) (sizeof(_fixed_trig_bitwise_rs_n_tab)
            / sizeof(short))
            && n >= _fixed_trig_bitwise_rs_n_tab[j + 1]; j++)
        ;
    return _fixed_trig_bitwise_rs_r_tab[j];
}

int
_fixed_tan_halfangle(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r)
{
    if (n < 1)
        return 0;

#if FLINT_BITS == 64
    /* fully specialized per-size implementations (trig_opt_<n>.c,
       emitted and tuned by dev/tune_fixed.py) */
    if (r == 0 && n <= 12)
    {
        static void (* const tab[])(nn_ptr, nn_ptr, nn_ptr, nn_srcptr) = {
            NULL,
            _fixed_trig_opt_1, _fixed_trig_opt_2, _fixed_trig_opt_3,
            _fixed_trig_opt_4, _fixed_trig_opt_5, _fixed_trig_opt_6,
            _fixed_trig_opt_7, _fixed_trig_opt_8, _fixed_trig_opt_9,
            _fixed_trig_opt_10, _fixed_trig_opt_11, _fixed_trig_opt_12
        };
        tab[n](ysin, ycos, ytan, x);
        return 1;
    }
#endif

    if (r == 0)
        r = fixed_trig_bitwise_rs_default_r(n);
    r = FLINT_MAX(r, 32);
    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

    /* no tabulated tangent series here: the residual contributes
       through sin and cos of t', whose quotient is never formed */
    _fixed_tan_halfangle_mid(ysin, ycos, ytan, x, n, r, NULL);
    return 1;
}

void
fixed_tan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)
{
    int ok;
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* n = 1 is unsupported on 32-bit limbs: the clamp in the
       half-angle path would drop r under the r >= 32 contract of the
       residual series */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

    /* the half-angle machinery covers every size on both word sizes
       (beyond the tabulated tangent series it takes sin and cos of the
       residual and never divides them: their ratio cancels) */
    ok = _fixed_tan_halfangle(NULL, NULL, res, x, n, r);
    FLINT_ASSERT(ok);
    (void) ok;
}
