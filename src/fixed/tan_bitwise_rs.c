
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

#include "tan_rotate.inc"

/* R = floor-or-slightly-below of (2^(128n-1) - 1) / S for S with the
   top bit set (S' = S 2^-64n in [1/2, 1)), i.e. R ~ 1/(2 S') scaled
   to n limbs, R < 2^(64n).  N holds the all-ones numerator and h is
   n limbs of scratch; both callers pass them.  Above the cutoff this
   runs fixed_inv_newton with one guard limb; its two-sided error
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
/* Evaluate only the sine series and recover the cosine through a
   squaring and a square root once the series gets long enough,
   mirroring EXP_USE_SINH: the sine series of the reduced argument
   t' < 2^-r drops one term per 2r bits, so FLINT_BITS n / (2 r)
   terms in all, and skipping the cosine pass saves one of the two
   rectangular-splitting sums against the fixed cost of the squaring
   and the square root.  50 terms matches the crossover measured on
   the reference machine (n = 150 with the default ladder, which
   puts n = 150 in the r = 192 band: 64 * 150 = 50 * 192); the
   development VM measured a much later crossover in context, in
   line with its inflated Newton cutoffs, and the exp analogue uses
   128 terms.  Retune with the Newton cutoffs. */
#ifndef FIXED_TRIG_SIN_SQRT_TERMS
#define FIXED_TRIG_SIN_SQRT_TERMS 50
#endif
#define TRIG_USE_SIN_SQRT(n, r) \
    (FLINT_BITS * (n) >= FIXED_TRIG_SIN_SQRT_TERMS * (slong) (r))

/* working precision from which that square root runs through
   fixed_sqrt_newton instead of mpn_sqrtrem (same shape as the exp
   sinh-path square root) */
#ifndef FIXED_TRIG_SQRT_NEWTON_CUTOFF
#define FIXED_TRIG_SQRT_NEWTON_CUTOFF 2000
#endif

/* precision from which the rotation product W is accumulated in the
   growing-extent form instead of one shift-and-add pair per factor
   over the full width */
#ifndef FIXED_TRIG_ROTATE_PROD_CUTOFF
#define FIXED_TRIG_ROTATE_PROD_CUTOFF 8
#endif

#ifndef FIXED_TRIG_RECIP_NEWTON_CUTOFF
#define FIXED_TRIG_RECIP_NEWTON_CUTOFF 40
#endif

static void
_fixed_recip_bn(nn_ptr R, nn_ptr h, nn_srcptr N, nn_srcptr S, slong n)
{
    if (n < FIXED_TRIG_RECIP_NEWTON_CUTOFF)
    {
        mpn_tdiv_qr(R, h, 0, N, 2 * n, S, n);
    }
    else
    {
        nn_ptr q;
        TMP_INIT;
        TMP_START;

        /* (q, n+3) ~ 1/S' with n+1 fraction limbs, 2 integral */
        q = TMP_ALLOC((n + 3) * sizeof(ulong));
        fixed_inv_newton(q, S, n, n + 1);

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

    nn_srcptr tab;
    tab = _fixed_atans_tab(n, r, &nc);

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
    num = _fixed_bitwise_reduce(t, n, r, 1, tab, nc, used);

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
    else if (n < FIXED_TRIG_ROTATE_PROD_CUTOFF)
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
    else
    {
        /* Growing-extent form of the same product, the trig analogue
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
        alloc = wn + 16;
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
        if (TRIG_USE_SIN_SQRT(n, r))
        {
            /* Evaluate only the sine series -- the analogue of the
               exp sinh path -- and recover

                   g = 1 - cos t' = 1 - sqrt(1 - sin^2 t')

               by a squaring and a square root, saving the cosine
               pass of fixed_sin_cos_rs (the powers of t'^2 are
               shared between the two sums, but each sum is a full
               rectangular-splitting pass).  The square root sits
               just below 1, so its derivative is ~1/2 and nothing
               is amplified: sqrhigh's few-ulp deficit on sin^2
               (which can only raise 1 - sin^2, lowering g), the
               square root's floor (raising g by at most one ulp)
               and the Newton route's ~2 compensated ulps leave g
               within a few 2^-64n ulps either way, well inside the
               budget the direct cosine's truncation already used. */
            nn_ptr u2 = N, rt = Q;      /* free until the divisions */

            fixed_sin_rs(ss, t, n);

            flint_mpn_sqrhigh(rt, ss, n);
            if (flint_mpn_zero_p(rt, n))
            {
                flint_mpn_zero(cc, n);  /* sin^2 below 2^-64n: g = 0 */
            }
            else if (n < FIXED_TRIG_SQRT_NEWTON_CUTOFF)
            {
                flint_mpn_zero(u2, n);
                mpn_neg(u2 + n, rt, n); /* (1 - sin^2) 2^(128n) */
                mpn_sqrtrem(rt, sc, u2, 2 * n);
                mpn_neg(cc, rt, n);     /* g 2^(64n) */
            }
            else
            {
                /* the Newton square root reads short input directly:
                   (1 - sin^2) as an n-limb fraction, no zero padding */
                mpn_neg(u2, rt, n);
                fixed_sqrt_newton(rt, u2, n, n + 2);
                if (rt[n + 2])
                    flint_mpn_zero(cc, n);  /* cos rounded to 1 */
                else
                    mpn_neg(cc, rt + 2, n);
            }
        }
        else
        {
            fixed_sin_cos_rs(ss, cc, t, n);

            if (cc[n])
                flint_mpn_zero(cc, n);      /* cos t' = 1: g = 0 */
            else
                mpn_neg(cc, cc, n);         /* g = 1 - cos t' */
        }

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
            _fixed_recip_bn(R, h, N, S, n);

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

            _fixed_recip_bn(R, h, N, W, n);

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
