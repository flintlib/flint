
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
#include "mp_real.h"
#include "impl.h"

/* R = floor-or-slightly-below of (2^(128n-1) - 1) / S for S with the
   top bit set (S' = S 2^-64n in [1/2, 1)), i.e. R ~ 1/(2 S') scaled
   to n limbs, R < 2^(64n).  N holds the all-ones numerator and h is
   n limbs of scratch; both callers pass them.  Above the cutoff this
   runs _mp_real_inv_newton with one guard limb; its two-sided error
   4 B^(-n-1) / S' <= 8 B^(-n-1), shifted into R's frame, stays below
   one R ulp, so truncating leaves R within 1 + eps ulps of
   floor(N / S) -- the callers' mulhigh truncation budgets absorb
   that like the tdiv floor.  When S' = 1/2 exactly the
   contract value 1/S' = 2 would overflow R; the compensated result
   clamps to all ones, matching the tdiv path.  The cutoff comes
   from in-context measurements (see dev/notes): in a tight loop the
   Newton inversion wins from a few dozen limbs, but embedded in the
   evaluation it loses 10-15% of the whole call at n = 512-1024 and
   only wins (3-6%) from about 1536 limbs, where the FFT range
   starts. */
#ifndef MP_REAL_TRIG_RECIP_NEWTON_CUTOFF
#define MP_REAL_TRIG_RECIP_NEWTON_CUTOFF 40
#endif

static void
_mp_real_recip_bn(nn_ptr R, nn_ptr h, nn_srcptr N, nn_srcptr S, slong n)
{
    if (n < MP_REAL_TRIG_RECIP_NEWTON_CUTOFF)
    {
        flint_mpn_tdiv_qr(R, h, N, 2 * n, S, n);
    }
    else
    {
        nn_ptr q;
        TMP_INIT;
        TMP_START;

        /* (q, n+3) ~ 1/S' with n+1 fraction limbs, 2 integral */
        q = TMP_ALLOC((n + 3) * sizeof(ulong));
        _mp_real_inv_newton(q, S, n, n + 1);

        /* R = (1/S') 2^(64n-1): shift q down 65 bits and truncate
           the guard limb; the two-sided error (< 8 half-ulps of q's
           bottom limb) leaves R within 1 + eps ulps of floor(N / S),
           which the callers' mulhigh truncation budgets absorb like
           the tdiv floor */
        mpn_rshift(q, q + 1, n + 2, 1);
        flint_mpn_copyi(R, q, n);

        /* q's remaining top limb is 1 iff 1/S' reached 2 (S' = 1/2):
           clamp to all ones as the all-ones numerator does */
        if (q[n] != 0)
            flint_mpn_store(R, n, ~UWORD(0));

        TMP_END;
    }
}

/* The tangent half-angle reconstruction, behind _mp_real_sin_cos_bitwise_rs
   and _mp_real_tan_bitwise_rs.

   The angle x/2 is reduced with the cached angles A_i = atan(2^-i),

       x/2 = sum_i d_i A_i + t',        0 <= t' < 2^-r,

   and the rotation W = prod (1 + d_i i 2^-i) turns the residual's
   (cos t', sin t') (or (1, tan t')) into (DE, TT), a positive real
   multiple of exp(i x/2): TT/DE = tan(x/2) = t.  Since TAN IS A RATIO the
   growth of |W| (below 1.1645) cancels and needs no normalization.  With

       A = TT DE,  B = DE^2,  C = TT^2,

   the half-angle formulas give sin x = 2A/(B + C), cos x = 1 - 2C/(B + C)
   and tan x = 2A/(B - C): one reciprocal of a divisor normalized to n
   limbs (its exponent folds into the final doubling shift) and one
   mulhigh per output, t itself never formed (_mp_real_tan_halfangle_finish).
   Note 0 <= t < tan(1/2) = 0.5464, so 2t may exceed one and the outputs
   carry a unit limb.

   Two reductions feed the finish: _mp_real_tan_halfangle_mid below, the
   greedy digits d_i in {0, 1} with limb arithmetic (32-bit machines, and
   explicit r below MP_REAL_TRIG_BITWISE_VEC_MIN limbs on 64-bit ones),
   and _mp_real_tan_halfangle_vec (bitwise_vec.c), signed digits on
   vectors.  With r = 0 at n <= MP_REAL_TRIG_OPT_R_MAX the per-size sine
   and cosine kernels serve instead (and tan their quotient,
   _mp_real_tan_sc). */

/* the outputs from TT = (T, n + 1), T[n] = 0, and DE = (D, n + 1),
   the imaginary and real parts of a positive multiple of
   exp(i x / 2), TT < 0.64 and DE in (0.72, 1.17); T and D are
   clobbered */
void
_mp_real_tan_halfangle_finish(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_ptr T, nn_ptr D, slong n)
{
    slong wn = n + 1;
    nn_ptr u, va, vb, N, Q, sc;
    TMP_INIT;

    TMP_START;
    u = TMP_ALLOC((n + n + n + (2 * n + 2) + (2 * n + 2) + (2 * n + 4))
        * sizeof(ulong));
    va = u + n;
    vb = va + n;
    N = vb + n;
    Q = N + (2 * n + 2);
    sc = Q + (2 * n + 2);

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
            _mp_real_recip_bn(R, h, N, S, n);

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

            _mp_real_recip_bn(R, h, N, W, n);

            /* tan = 2A/(B - C) = 2^(2+e) A R < 1.56 */
            flint_mpn_mulhigh_n(ytan, A, R, n);
            ytan[n] = 0;
            mpn_lshift(ytan, ytan, n + 1, 2 + e);
        }
    }

    TMP_END;
}

void
_mp_real_tan_halfangle_mid(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r)
{
    slong wn = n + 1, nc, num, j;
    nn_ptr t, wx, wy, va, vb, T, D, ss, cc;
    slong * used;
    TMP_INIT;

    nn_srcptr tab;
    tab = _mp_real_atans_tab(n, r, &nc);

    TMP_START;
    t = TMP_ALLOC((n + wn + wn + (wn + 1) + n + wn + wn
        + (n + 2) + (n + 2)) * sizeof(ulong));
    wx = t + n;
    wy = wx + wn;
    va = wy + wn;               /* wn + 1: the expansion of alpha */
    vb = va + wn + 1;
    T = vb + n;                 /* wn: TT */
    D = T + wn;                 /* wn: DE */
    ss = D + wn;                /* n+2 each: sin, cos of the residual */
    cc = ss + (n + 2);
    used = TMP_ALLOC(MP_REAL_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    /* reduce x/2; the residual t' is wanted as it stands, so unlike the
       sin/cos reduction there is no doubling and no extra index */
    mpn_rshift(t, x, n, 1);
    num = _mp_real_bitwise_reduce(t, n, r, 1, tab, nc, used);

    /* W = prod (1 + i 2^-i) */
    {
        /* Growing-extent form of the rotation product, the trig analogue
           of the exp product reconstruction: with the used indices
           ascending write W = (1 - alpha) + i b, alpha, b >= 0
           (alpha starts at 0 and only grows; b stays below 1 since
           the total angle is under one radian), both lsb-anchored as
           integers at a common extent B that deepens by i per factor:

               alpha' = (alpha << i) + b
               b'     = (b << i) + 2^B - alpha

           using the old values on the right.  The extent grows only
           quadratically (about the sum of the used indices), so
           several factors fit in one limb, then two, three, ...;
           limbs deeper than 64 (wn + 1) bits are dropped as the
           buffers grow, each drop losing under one bit two guard
           limbs down.  Unlike exp there is no final multiplication
           at all: wx = 1 - alpha and wy = b are expanded into their
           dense frames in one pass each. */
        slong capd = FLINT_BITS * (wn + 1);
        slong B = 0, al = 1, bl = 1, alloc;
        nn_ptr aa, at, ba, bt;
        TMP_INIT;

        TMP_START;
        /* extent below capd + 63 + max_index bits; a fixed "+ 16"
           slack overflowed on 32-bit limbs at large r, exactly as
           in the exp product reconstruction */
        alloc = wn + 4 + ((num > 0) ? used[num - 1] : 0) / FLINT_BITS;
        aa = TMP_ALLOC(4 * alloc * sizeof(ulong));
        at = aa + alloc;
        ba = at + alloc;
        bt = ba + alloc;
        aa[0] = 0;
        ba[0] = 0;

        for (j = 0; j < num; j++)
        {
            slong ii = used[j], q = ii / FLINT_BITS;
            int sb = (int) (ii - q * FLINT_BITS);
            slong la, lb, d;
            ulong cy;

            /* at = (aa << ii) + ba */
            la = al + q + (sb != 0);
            flint_mpn_zero(at, q);
            if (sb)
                at[al + q] = mpn_lshift(at + q, aa, al, sb);
            else
                flint_mpn_copyi(at + q, aa, al);
            while (la > 1 && at[la - 1] == 0)
                la--;
            if (la >= bl)
            {
                cy = mpn_add(at, at, la, ba, bl);
                at[la] = cy;
                la += (cy != 0);
            }
            else
            {
                cy = mpn_add(at, ba, bl, at, la);
                at[bl] = cy;
                la = bl + (cy != 0);
            }

            /* bt = (ba << ii) + 2^B - aa */
            lb = bl + q + (sb != 0);
            flint_mpn_zero(bt, q);
            if (sb)
                bt[bl + q] = mpn_lshift(bt + q, ba, bl, sb);
            else
                flint_mpn_copyi(bt + q, ba, bl);
            while (lb > 1 && bt[lb - 1] == 0)
                lb--;
            {
                slong pb = B / FLINT_BITS;
                if (pb >= lb)
                {
                    flint_mpn_zero(bt + lb, pb - lb + 1);
                    lb = pb + 1;
                }
                cy = mpn_add_1(bt + pb, bt + pb, lb - pb,
                    UWORD(1) << (B % FLINT_BITS));
                bt[lb] = cy;
                lb += (cy != 0);
            }
            if (al > 0 && !(al == 1 && aa[0] == 0))
                mpn_sub(bt, bt, lb, aa, al);
            while (lb > 1 && bt[lb - 1] == 0)
                lb--;

            { nn_ptr u_; u_ = aa; aa = at; at = u_;
                         u_ = ba; ba = bt; bt = u_; }
            al = la;
            bl = lb;
            B += ii;

            d = (B - capd) / FLINT_BITS;
            if (d > 0)
            {
                if (d < al) { flint_mpn_copyi(aa, aa + d, al - d); al -= d; }
                else { aa[0] = 0; al = 1; }
                if (d < bl) { flint_mpn_copyi(ba, ba + d, bl - d); bl -= d; }
                else { ba[0] = 0; bl = 1; }
                B -= FLINT_BITS * d;
            }
        }

        /* expand: wy = b 2^(64 n), wx = 2^(64 n) - alpha 2^(64 n) */
        {
            slong sh = FLINT_BITS * n - B;

            flint_mpn_zero(wy, wn);
            flint_mpn_zero(va, wn + 1);
            if (sh >= 0)
            {
                slong q = sh / FLINT_BITS;
                int b2 = (int) (sh % FLINT_BITS);
                if (b2)
                {
                    if (bl + q < wn + 1)
                        wy[bl + q] = mpn_lshift(wy + q, ba, bl, b2);
                    else
                        mpn_lshift(wy + q, ba, bl, b2);
                    va[al + q] = mpn_lshift(va + q, aa, al, b2);
                }
                else
                {
                    flint_mpn_copyi(wy + q, ba, bl);
                    flint_mpn_copyi(va + q, aa, al);
                }
            }
            else
            {
                slong q = (-sh) / FLINT_BITS;
                int b2 = (int) ((-sh) % FLINT_BITS);
                if (b2)
                {
                    if (bl - q > 0)
                        mpn_rshift(wy, ba + q, bl - q, b2);
                    if (al - q > 0)
                        mpn_rshift(va, aa + q, al - q, b2);
                }
                else
                {
                    if (bl - q > 0)
                        flint_mpn_copyi(wy, ba + q, bl - q);
                    if (al - q > 0)
                        flint_mpn_copyi(va, aa + q, al - q);
                }
            }
            flint_mpn_zero(wx, n);
            wx[n] = 1;
            mpn_sub(wx, wx, wn, va, wn);
        }
        TMP_END;
    }

    /* TT and DE, the imaginary and real parts of W (cos t' + i sin t')
       up to positive real factors that CANCEL in every ratio below just
       as |W| does: TT/DE is tan(x/2).  TT, at most
       tan(1/2) = 0.5464 times DE in (0.72, 1.17), stays below 0.64 and
       never uses its unit limb.  With cos t' = 1 - g, g < 2^(-2r-1),

           TT = wy cos t' + wx sin t' = wy + wx ss - wy g,
           DE = wx cos t' - wy sin t' = wx - wx g - wy ss,

       four mulhighs against the tiny ss and g; no difference can go
       negative: TT's subtrahend wy g is below wy, DE's two are below
       wx. */
    {
        /* sin t' and g = 1 - cos t' of the residual
           (_mp_real_sin_cos_reduced also owns the tuned series-vs-sqrt
           choice and the bit-burst modes for large n / r); r >= 32 here
           per _mp_real_tan_halfangle */
        _mp_real_sin_cos_reduced(ss, cc, NULL, t, n, (flint_bitcnt_t) r, 0);

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

    _mp_real_tan_halfangle_finish(ysin, ycos, ytan, T, D, n);

    TMP_END;
}

/* smallest n at which each reduction parameter becomes optimal for
   the shared half-angle path of _mp_real_sin_cos_bitwise_rs and
   _mp_real_tan_bitwise_rs beyond the range of the tapered tangent
   series: interleaved timings of every ladder value at 17 sizes from
   84 to 1300 limbs (x86-64; the optimum flat to a few percent over
   one or two ladder steps).  The earlier table from
   src/mp_real/tune/tune-bitwise-r.c chose r one or two steps too
   large in places (2-9% slower). */
static const int _mp_real_trig_bitwise_rs_r_tab[] =
    {32, 64, 128, 192, 256, 320, 384, 448, 512};
static const short _mp_real_trig_bitwise_rs_n_tab[] =
    {1, 11, 40, 130, 330, 520, 700, 950, 1600};

#if FLINT_BITS == 64
/* the same for _mp_real_tan_halfangle_vec, from
   MP_REAL_TRIG_BITWISE_VEC_MIN limbs (its own sweep, x86-64 AVX2:
   every ladder value at 27 sizes from 4 to 600 limbs, the minimum of
   three interleaved runs, the optimum flat to a few percent over one or
   two steps; of sin and cos, which share it with tan) */
static const int _mp_real_trig_bitwise_vec_r_tab[] =
    {48, 64, 96, 128, 192, 256, 320, 384, 448, 512, 640};
static const short _mp_real_trig_bitwise_vec_n_tab[] =
    {1, 28, 32, 44, 90, 120, 150, 230, 280, 380, 420};
#endif

int
_mp_real_trig_bitwise_rs_default_r(slong n)
{
    slong j;

#if FLINT_BITS == 64
    /* the reduction parameters of the per-size kernels that r = 0
       selects at these sizes (see impl.h) */
    static const int opt_r[] = MP_REAL_TRIG_OPT_R;
    if (n <= MP_REAL_TRIG_OPT_R_MAX)
        return opt_r[n];

    if (n >= MP_REAL_TRIG_BITWISE_VEC_MIN)
    {
        for (j = 0; j + 1 < (slong) (sizeof(_mp_real_trig_bitwise_vec_n_tab)
                / sizeof(short))
                && n >= _mp_real_trig_bitwise_vec_n_tab[j + 1]; j++)
            ;
        return _mp_real_trig_bitwise_vec_r_tab[j];
    }
#endif

    for (j = 0; j + 1 < (slong) (sizeof(_mp_real_trig_bitwise_rs_n_tab)
            / sizeof(short))
            && n >= _mp_real_trig_bitwise_rs_n_tab[j + 1]; j++)
        ;
    return _mp_real_trig_bitwise_rs_r_tab[j];
}

/* (y, n + 1) = tan v from the sine and cosine kernels and one division:
   with the sine and cosine within E ulps and cos v > cos 1 > 0.5403, the
   quotient is within (1 + tan 1) E / 0.5403 + 1 < 4.74 E + 1 ulps */
#if FLINT_BITS == 64
void
_mp_real_tan_sc(nn_ptr y, ulong * err, nn_srcptr v, slong n)
{
    nn_ptr s, c;
    ulong e;
    TMP_INIT;

    TMP_START;
    s = TMP_ALLOC((n + 1 + n + 1) * sizeof(ulong));
    c = s + n + 1;

    /* the kernels cover every n of the contract */
    if (!_mp_real_sin_cos_opt(s, c, &e, v, n))
        flint_throw(FLINT_ERROR, "_mp_real_tan_sc: n = %wd beyond the "
            "sine and cosine kernels\n", n);

    if (c[n] != 0)
    {
        /* cos v = 1 to the working precision */
        flint_mpn_copyi(y, s, n);
        y[n] = 0;
    }
    else if (n == 1)
    {
        /* one hardware division, the units limb from s >= c */
        ulong s0 = s[0], c0 = c[0], r;
        y[1] = (s0 >= c0);
        if (s0 >= c0)
            s0 -= c0;
        udiv_qrnnd(y[0], r, s0, UWORD(0), c0);
        (void) r;
    }
    else
    {
        /* c >= 0.54: normalized at n limbs; the quotient takes n + 1,
           within one ulp of the floor either way */
        flint_mpn_divapprox_fraction(y, s, n, c, n, n);
    }

    if (err != NULL)
        *err = (19 * e) / 4 + 2;

    TMP_END;
}
#endif

int
_mp_real_tan_halfangle(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r)
{
    if (n < 1)
        return 0;

    if (r == 0)
        r = _mp_real_trig_bitwise_rs_default_r(n);
    r = FLINT_MAX(r, 32);
    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

#if FLINT_BITS == 64
    if (n >= MP_REAL_TRIG_BITWISE_VEC_MIN)
    {
        _mp_real_tan_halfangle_vec(ysin, ycos, ytan, x, n, r);
        return 1;
    }
#endif

    _mp_real_tan_halfangle_mid(ysin, ycos, ytan, x, n, r);
    return 1;
}

FLINT_FORCE_INLINE void
_mp_real_tan_bitwise_rs_noerr(nn_ptr res, nn_srcptr x, slong n, int r)
{
    int ok;
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* n = 1 is unsupported on 32-bit limbs: the clamp in the
       half-angle path would drop r under the r >= 32 contract of the
       residual series */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

    ok = _mp_real_tan_halfangle(NULL, NULL, res, x, n, r);
    FLINT_ASSERT(ok);
    (void) ok;
}

void
_mp_real_tan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)
{
#if FLINT_BITS == 64
    /* r = 0 at the sizes of the per-size sine and cosine kernels: those
       and one division */
    if (r == 0 && n <= MP_REAL_TRIG_OPT_R_MAX)
    {
        _mp_real_tan_sc(res, err, x, n);
        return;
    }
#endif

    _mp_real_tan_bitwise_rs_noerr(res, x, n, r);
    if (err != NULL)
        *err = 8 * (r ? (ulong) r : (ulong) _mp_real_trig_bitwise_rs_default_r(n)) + 256;
}
