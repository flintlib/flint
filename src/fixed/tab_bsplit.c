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
#include "fixed.h"

/* Binary splitting for the table values log(1 + 2^-i) and
   atan(2^-i), in pure mpn arithmetic.

   Both series have dyadic ratios (2^-i resp. 2^-2i), except that the
   logarithm switches to log(1 + 2^-i) = 2 atanh(1/q), q = 2^(i+1)+1,
   below the direct-series cutoff; there the ratio q^-2 is not dyadic,
   but multiplication by q is a shift and an addition.

   INTEGERS AS LIMB-RADIX FLOATS.  The split products are carried as
   normalized triples (d, size, exp): value = (d, size) 2^(64 exp),
   top limb nonzero, trailing zero LIMBS stripped into exp (bit-level
   trailing zeros stay in the mantissa).  Whenever a mantissa would
   exceed the cap (the target size plus a few limbs of margin), it is
   truncated and the dropped limbs move into the exponent.

   ONE-SIDEDNESS.  Every truncation on the numerator side rounds
   DOWN and every one on the denominator side rounds UP, the final
   division truncates, and the omitted series tail is positive (the
   term count is even-aligned, so for the alternating series the
   first omitted term has a plus sign).  The result is therefore at
   most the true value, short by a couple of ulps -- exactly the
   contract of the reduction tables.

   SIGNS.  The alternating series are evaluated over blocks aligned
   to even indices: every block sum and every merge is then a sum of
   positive quantities, and only the iterative basecase subtracts
   (with a positive running value, since the terms decrease).

   BASECASE.  Instead of recursing to length 1 or 2, blocks of
   TAB_BS_BLK terms are accumulated iteratively with linear
   operations: with V(b) = 0 and

       V(j) = 1/d_j (+/-) x V(j+1),

   maintaining V = N / (D 2^F) (dyadic x = 2^-e) gives

       N <- (D << (F + e - v_j)) (-/+) o_j N,      D <- o_j D,

   where d_j = o_j 2^(v_j) with o_j odd -- a shift of the small
   denominator product plus an mpn_(add/sub)mul_1.  The atanh
   basecase is the same with the shift of D replaced by
   multiplication by q^2 = (D << (2i+2)) + (D << (i+2)) + D. */

#define TAB_BS_BLK 8

/* dyadic_finish divides by Newton iteration (fixed_div_newton)
   instead of mpn_tdiv_qr when the denominator has at least this many
   limbs AND at least a third of the numerator's: with a short
   divisor GMP's division exploits the shape while the Newton
   division always pays for a full-precision reciprocal (measured
   in situ: 1.15-1.85x faster above the gate, as low as 0.5x below) */
#ifndef TAB_BS_DIV_NEWTON_CUTOFF
#define TAB_BS_DIV_NEWTON_CUTOFF 40
#endif

/* extra buffer limbs beyond the truncation cap: a basecase block
   spans up to TAB_BS_BLK terms of e bits each, plus small change */
#define BUFSLACK(e) ((TAB_BS_BLK * (e)) / FLINT_BITS + 8)

/* triple (d, n, e): value = (d, n) 2^(FLINT_BITS e); n = 0 means 0 */
typedef struct
{
    nn_ptr d;
    slong n;
    slong e;
}
fx_t;

/* strip leading zero limbs and move trailing zero limbs into e */
static void
fx_norm(fx_t * x)
{
    slong t;

    while (x->n > 0 && x->d[x->n - 1] == 0)
        x->n--;

    if (x->n == 0)
    {
        x->e = 0;
        return;
    }

    for (t = 0; x->d[t] == 0; t++)
        ;

    if (t > 0)
    {
        flint_mpn_copyi(x->d, x->d + t, x->n - t);
        x->n -= t;
        x->e += t;
    }
}

/* truncate the mantissa to at most cap limbs; round_up adds one ulp
   (of the truncated mantissa) when any dropped limb was nonzero */
static void
fx_trunc(fx_t * x, slong cap, int round_up)
{
    slong drop = x->n - cap;
    int inexact = 0;
    slong t;

    if (drop <= 0)
        return;

    if (round_up)
        for (t = 0; t < drop && !inexact; t++)
            inexact = (x->d[t] != 0);

    flint_mpn_copyi(x->d, x->d + drop, cap);
    x->n = cap;
    x->e += drop;

    if (inexact)
    {
        if (mpn_add_1(x->d, x->d, cap, 1))
        {
            /* carry out: value 2^(64 cap) */
            x->d[0] = 1;
            x->n = 1;
            x->e += cap;
            return;
        }
    }

    fx_norm(x);
}

/* dst = a * b, truncated to cap; dst->d must hold a->n + b->n limbs
   (callers reserve cap + 4 and keep operands at most cap + 1 after
   their own truncation, with muls going through scratch when both
   operands are large) */
static void
fx_mul(fx_t * dst, const fx_t * a, const fx_t * b, slong cap,
    int round_up)
{
    slong ps = a->n + b->n;

    if (a->n == 0 || b->n == 0)
    {
        dst->n = 0;
        dst->e = 0;
        return;
    }

    /* When the full product exceeds the cap and the operands are
       near-balanced -- the steady state at the top of the splitting
       tree, where truncation pins both at cap -- compute only the
       kept window with the windowed middle product instead of
       forming and discarding the low limbs.  flint_mpn_mulmid is a
       LOWER approximation: it never exceeds the exact window, and
       the deficit is bounded by min(an, bn, zlo) 2^64 window units
       -- up to ~min ulps of the limb ABOVE the window bottom, a
       bound that all-ones operands nearly attain -- NOT by min + 1
       ulps of the bottom limb as an earlier version of this comment
       claimed.  The round-down case still needs no correction (the
       deficit points down and sits two limbs below the entry
       guard), while the round-up case must add min + 2 at limb 1,
       covering the full deficit; the overshoot, of the same
       magnitude, lands equally far below the guard.  Measured
       1.1-1.4x faster than
       flint_mpn_mul at these shapes from 64 to 10240 limbs; the
       skewed shapes stay on the full product, where the mulmid
       dispatcher is not yet competitive below ~4096 limbs. */
    if (ps > cap && FLINT_MIN(a->n, b->n) >= cap - 4)
    {
        flint_mpn_mulmid(dst->d, a->d, a->n, b->d, b->n,
            ps - cap, ps);
        dst->n = cap;
        dst->e = a->e + b->e + (ps - cap);

        if (round_up)
        {
            if (mpn_add_1(dst->d + 1, dst->d + 1, cap - 1,
                (ulong) FLINT_MIN(a->n, b->n) + 2))
            {
                dst->d[0] = 1;
                dst->n = 1;
                dst->e += cap;
                return;
            }
        }

        fx_norm(dst);
        return;
    }

    if (a->n >= b->n)
        flint_mpn_mul(dst->d, a->d, a->n, b->d, b->n);
    else
        flint_mpn_mul(dst->d, b->d, b->n, a->d, a->n);

    dst->n = a->n + b->n;
    dst->e = a->e + b->e;
    fx_norm(dst);
    fx_trunc(dst, cap, round_up);
}

/* x <<= s bits, in place; d must hold n + s/64 + 1 limbs */
static void
fx_shift_up(fx_t * x, slong s)
{
    slong q = s / FLINT_BITS;
    int b = (int) (s % FLINT_BITS);

    if (x->n == 0 || s == 0)
    {
        x->e += (x->n == 0) ? 0 : q;
        if (x->n != 0 && b != 0)
        {
            ulong cy = mpn_lshift(x->d, x->d, x->n, b);
            if (cy)
                x->d[x->n++] = cy;
        }
        return;
    }

    x->e += q;
    if (b != 0)
    {
        ulong cy = mpn_lshift(x->d, x->d, x->n, b);
        if (cy)
            x->d[x->n++] = cy;
        fx_norm(x);
    }
}

/* dst = a + b (both nonnegative), truncated to cap rounding down;
   dst->d must hold cap + 3 limbs; dst may alias neither operand */
static void
fx_add(fx_t * dst, const fx_t * a, const fx_t * b, slong cap)
{
    const fx_t * lo, * hi;
    slong off, top, bot, span;

    if (a->n == 0) { flint_mpn_copyi(dst->d, b->d, b->n);
                     dst->n = b->n; dst->e = b->e; return; }
    if (b->n == 0) { flint_mpn_copyi(dst->d, a->d, a->n);
                     dst->n = a->n; dst->e = a->e; return; }

    if (a->e <= b->e) { lo = a; hi = b; }
    else              { lo = b; hi = a; }

    top = FLINT_MAX(lo->e + lo->n, hi->e + hi->n);
    bot = lo->e;
    span = top - bot;

    if (span > cap + 2)
    {
        /* the low limbs cannot survive truncation: drop them
           (rounding down) and re-anchor */
        bot = top - (cap + 2);
        span = cap + 2;
    }

    /* assemble hi at its offset, over [bot, top) */
    off = hi->e - bot;
    if (off >= 0)
    {
        flint_mpn_zero(dst->d, off);
        flint_mpn_copyi(dst->d + off, hi->d, hi->n);
    }
    else
    {
        /* hi partially below the window (possible only when lo was
           re-anchored above hi->e; then hi == the smaller-exponent...
           cannot happen: bot <= lo->e <= hi->e) */
        FLINT_ASSERT(0);
    }
    flint_mpn_zero(dst->d + off + hi->n, span - off - hi->n);
    dst->n = span;
    dst->e = bot;

    /* add lo where it overlaps the window */
    off = lo->e - bot;
    if (off + lo->n > 0)
    {
        slong skip = (off < 0) ? -off : 0;    /* skip == 0 here */
        ulong cy = mpn_add(dst->d + off + skip, dst->d + off + skip,
            dst->n - off - skip, lo->d + skip, lo->n - skip);
        if (cy)
        {
            dst->d[dst->n++] = cy;            /* room: cap + 3 */
        }
    }

    fx_norm(dst);
    fx_trunc(dst, cap, 0);
}

/* ---------------------------------------------------------------- */
/* dyadic engine: S(a, b) = sum_{j=a}^{b-1} (-1)^(j-a) x^(j-a) / d_j
   with x = 2^-e and d_j = odd(base + step j) 2^(v_j); a and b even.
   Output: S = N / (D 2^F), N and D nonnegative, D odd, F in BITS.
   N and D are written into caller-provided fx buffers of cap + 4
   limbs; scratch provides working space for this frame and every
   frame below it. */

typedef struct
{
    slong e;             /* bit shift per term */
    ulong dbase, dstep;  /* d_j = dbase + dstep * j */
    slong cap;
}
dyadic_t;

static void
dyadic_bs(fx_t * N, fx_t * D, slong * F, slong a, slong b,
    const dyadic_t * par, nn_ptr scratch, slong scratch_alloc)
{
    slong cap = par->cap;

    if (b - a <= TAB_BS_BLK)
    {
        /* iterative basecase, top down */
        slong j;
        ulong dj, oj;
        int vj;
        fx_t T;

        T.d = scratch;              /* cap + BUFSLACK(e) limbs */

        j = b - 1;
        dj = par->dbase + par->dstep * (ulong) j;
        oj = dj >> flint_ctz(dj);
        vj = flint_ctz(dj);

        N->d[0] = 1; N->n = 1; N->e = 0;
        D->d[0] = oj; D->n = 1; D->e = 0;
        *F = vj;

        for (j = b - 2; j >= a; j--)
        {
            slong s;

            dj = par->dbase + par->dstep * (ulong) j;
            vj = flint_ctz(dj);
            oj = dj >> vj;

            /* T = D << (F + e - v_j)  (in bits; exponent limb-split) */
            s = *F + par->e - vj;
            FLINT_ASSERT(s >= 0);
            flint_mpn_copyi(T.d, D->d, D->n);
            T.n = D->n; T.e = D->e;
            fx_shift_up(&T, s);

            /* N = T - o_j N  (positive: alternating decreasing) */
            {
                /* align: N scaled by o_j, then subtract from T at
                   exponent offsets.  Materialize both at T's window:
                   N has smaller magnitude. */
                slong off = N->e - T.e;
                ulong cy;

                /* multiply N by o_j in place (room: cap + 4) */
                cy = mpn_mul_1(N->d, N->d, N->n, oj);
                if (cy) N->d[N->n++] = cy;

                if (off < 0)
                {
                    /* extend T downward so the subtraction is exact
                       (the buffers carry enough slack: a block spans
                       at most BLK terms of e bits each) */
                    slong ext = -off;
                    slong t;
                    for (t = T.n - 1; t >= 0; t--)
                        T.d[t + ext] = T.d[t];
                    flint_mpn_zero(T.d, ext);
                    T.n += ext;
                    T.e -= ext;
                    off = 0;
                }

                FLINT_ASSERT(off + N->n <= T.n);
                cy = mpn_sub(T.d + off, T.d + off, T.n - off,
                    N->d, N->n);
                FLINT_ASSERT(cy == 0);
            }

            flint_mpn_copyi(N->d, T.d, T.n);
            N->n = T.n; N->e = T.e;
            fx_norm(N);
            fx_trunc(N, cap, 0);

            /* D <- o_j D */
            {
                ulong cy = mpn_mul_1(D->d, D->d, D->n, oj);
                if (cy) D->d[D->n++] = cy;
                fx_trunc(D, cap, 1);
            }

            *F = *F + par->e;
        }
        return;
    }

    {
        slong m = a + (((b - a) / 2 + 1) & ~(slong) 1);
        fx_t N2, D2, T1, T2;
        slong F2, s;
        nn_ptr sc = scratch;

        slong bufn = cap + BUFSLACK(par->e);

        N2.d = sc; sc += bufn;
        D2.d = sc; sc += bufn;
        T1.d = sc; sc += 2 * bufn;
        T2.d = sc; sc += 2 * bufn;

        dyadic_bs(N, D, F, a, m, par, sc,
            scratch_alloc - (sc - scratch));
        dyadic_bs(&N2, &D2, &F2, m, b, par, sc,
            scratch_alloc - (sc - scratch));

        /* S = N1/(D1 2^F1) + 2^(-e (m-a)) N2/(D2 2^F2)
             = [N1 D2 2^(F2 + e(m-a) - F1) + N2 D1]
               / (D1 D2 2^(F2 + e(m-a))) */
        fx_mul(&T1, N, &D2, cap, 0);
        fx_mul(&T2, &N2, D, cap, 0);

        s = F2 + par->e * (m - a) - *F;
        if (s >= 0)
            fx_shift_up(&T1, s);
        else
            fx_shift_up(&T2, -s);             /* rebalance: multiply
                the OTHER side; denominator exponent adjusts below */

        fx_add(N, &T1, &T2, cap);
        fx_mul(&T1, D, &D2, cap, 1);
        flint_mpn_copyi(D->d, T1.d, T1.n);
        D->n = T1.n; D->e = T1.e;

        *F = F2 + par->e * (m - a) + ((s < 0) ? -s : 0);
    }
}

/* result = floor-ish of (N / (D 2^F)) * 2^(-outer) * 2^(64 n),
   one-sided (at most the true value); res gets n limbs */
static void
dyadic_finish(nn_ptr res, slong n, const fx_t * N, const fx_t * D,
    slong F, slong outer)
{
    slong s, num_n, den_n, qn;
    nn_ptr num, den, q, r;
    TMP_INIT;

    /* target = N 2^(64 N->e) 2^(64 n - F - outer) / (D 2^(64 D->e)) */
    s = FLINT_BITS * (n + N->e - D->e) - F - outer;

    TMP_START;

    if (s >= 0)
    {
        num_n = N->n + s / FLINT_BITS + 2;
        num = TMP_ALLOC(num_n * sizeof(ulong));
        flint_mpn_zero(num, num_n);
        flint_mpn_copyi(num + s / FLINT_BITS, N->d, N->n);
        if (s % FLINT_BITS)
            mpn_lshift(num + s / FLINT_BITS, num + s / FLINT_BITS,
                N->n + 1, (int) (s % FLINT_BITS));
        den_n = D->n;
        den = (nn_ptr) D->d;
    }
    else
    {
        num_n = N->n;
        num = (nn_ptr) N->d;
        den_n = D->n + (-s) / FLINT_BITS + 2;
        den = TMP_ALLOC(den_n * sizeof(ulong));
        flint_mpn_zero(den, den_n);
        flint_mpn_copyi(den + (-s) / FLINT_BITS, D->d, D->n);
        if ((-s) % FLINT_BITS)
            mpn_lshift(den + (-s) / FLINT_BITS,
                den + (-s) / FLINT_BITS, D->n + 1,
                (int) ((-s) % FLINT_BITS));
    }

    while (num_n > 0 && num[num_n - 1] == 0) num_n--;
    while (den_n > 0 && den[den_n - 1] == 0) den_n--;

    flint_mpn_zero(res, n);

    if (num_n >= den_n && num_n > 0)
    {
        if (den_n >= TAB_BS_DIV_NEWTON_CUTOFF && 3 * den_n >= num_n)
        {
            /* Newton division instead of mpn_tdiv_qr: view both
               operands as fractions, so the wanted quotient is
               (num/B^num_n) / (den/B^den_n) * B^k with
               k = num_n - den_n, and compute it with k + 3 fraction
               limbs; the numerator's long run of low zero limbs
               (from the shift) is skipped, which leaves the value
               unchanged.  fixed_div_newton's error is two-sided, at
               most 4 B^(3-k) / (den/B^den_n) <= 4 B^(2-k) -- far
               below B^-1 of an output ulp -- so truncating the three
               guard limbs yields floor(num/den) or one below it, and
               the entry stays one-sided like the tdiv floor. */
            slong k = num_n - den_n;
            slong nd = k + 3;
            slong j, tz = 0;

            while (tz < num_n && num[tz] == 0)
                tz++;

            q = TMP_ALLOC((nd + 2) * sizeof(ulong));
            fixed_div_newton(q, num + tz, num_n - tz, den, den_n, nd);

            for (j = 0; j < n; j++)
                res[j] = (3 + j < nd + 2) ? q[3 + j] : 0;
        }
        else
        {
            qn = num_n - den_n + 1;
            q = TMP_ALLOC((qn + 1) * sizeof(ulong));
            r = TMP_ALLOC(den_n * sizeof(ulong));
            mpn_tdiv_qr(q, r, 0, num, num_n, den, den_n);
            while (qn > 0 && q[qn - 1] == 0) qn--;
            FLINT_ASSERT(qn <= n);
            flint_mpn_copyi(res, q, qn);
        }
    }

    TMP_END;
}

void
fixed_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
{
    dyadic_t par;
    slong Nterms, F, prec, depth, alloc;
    fx_t N, D;
    nn_ptr buf, scratch;

    FLINT_ASSERT(i >= 1);
    FLINT_ASSERT(n >= 1);

    prec = FLINT_BITS * n + 8;
    Nterms = prec / (2 * (slong) i) + 2;
    Nterms = (Nterms + 1) & ~(slong) 1;       /* even */

    par.e = 2 * (slong) i;
    par.dbase = 1;
    par.dstep = 2;                            /* d_j = 2j + 1 */
    par.cap = n + 3;

    for (depth = 1; (TAB_BS_BLK << depth) < Nterms; depth++)
        ;
    /* every buffer holds a full basecase block span (BLK terms of
       e bits) on top of the cap */
    alloc = (depth + 2) * (6 * (par.cap + BUFSLACK(par.e)) + 40)
        + 4 * (par.cap + BUFSLACK(par.e)) + 40;

    buf = flint_malloc((2 * (par.cap + BUFSLACK(par.e)) + alloc)
        * sizeof(ulong));
    N.d = buf;
    D.d = buf + (par.cap + BUFSLACK(par.e));
    scratch = buf + 2 * (par.cap + BUFSLACK(par.e));

    dyadic_bs(&N, &D, &F, 0, Nterms, &par, scratch, alloc);
    dyadic_finish(res, n, &N, &D, F, (slong) i);

    flint_free(buf);
}

/* ---------------------------------------------------------------- */
/* atanh engine, computing
   W(a, b) = q^-2 sum_{j=a}^{b-1} q^(-2(j-a)) / (2j+1),
   q = 2^(i+1) + 1, all terms positive.  Output W = N / DEN with the
   denominator fully materialized as DEN = prod d_j q^(2(b-a)) --
   one q^2 per TERM, not per term less one -- plus the bare odd part
   D = prod d_j.  With that convention ranges compose without any
   q-power at all:

       W(a,b) = W(a,m) + q^(-2(m-a)) W(m,b)
              = (N1 DEN2 + N2 D1) / (DEN1 DEN2),

   four products per merge.  The q^2 factors enter only at the
   leaves, one linear multiplication
   x q^2 = (x << (2i+2)) + (x << (i+2)) + x per term (a single
   mpn_mul_1 whenever q^2 fits in a limb, i.e. 2i + 2 < FLINT_BITS);
   no full-size q-power is ever formed.

   Rounding: N, D and DEN all truncate DOWN (DEN appears in merge
   numerators, so it cannot round up); the final denominator is
   instead INFLATED by a rigorous bound on its rounded-down history,
   which restores one-sidedness at a relative cost of
   2^(-64 (cap-1) + ~13), far below an output ulp. */

/* x <- x q^2 = (x << (2i+2)) + (x << (i+2)) + x; scratch: room for
   x->n + (2i+2)/64 + 4 limbs plus a shifted copy of x */
static void
fx_mul_q2(fx_t * x, ulong i, slong cap, int round_up, nn_ptr scratch)
{
    fx_t T;
    slong sh = 2 * (slong) i + 2, q, off;
    int b;
    ulong cy;

    if (x->n == 0)
        return;

    /* q^2 fits in a limb whenever 2i + 2 < FLINT_BITS -- in
       particular for every i below the direct-series cutoff on
       64-bit -- and the three shift-and-adds collapse to a single
       mpn_mul_1 pass.  (The product of odd values stays odd, so no
       trailing limb can appear and fx_norm is unnecessary.) */
    if (2 * i + 2 < FLINT_BITS)
    {
        ulong qsq = (UWORD(1) << (2 * i + 2))
            + (UWORD(1) << (i + 2)) + 1;
        cy = mpn_mul_1(x->d, x->d, x->n, qsq);
        if (cy)
            x->d[x->n++] = cy;
        fx_trunc(x, cap, round_up);
        return;
    }

    T.d = scratch;
    flint_mpn_copyi(T.d, x->d, x->n);
    T.n = x->n; T.e = x->e;
    fx_shift_up(&T, sh);

    /* T += x << (i+2): limb offset relative to T.e */
    sh = (slong) i + 2;
    q = sh / FLINT_BITS + x->e;
    b = (int) (sh % FLINT_BITS);
    off = q - T.e;
    if (off < 0)
    {
        slong ext = -off, t;
        for (t = T.n - 1; t >= 0; t--)
            T.d[t + ext] = T.d[t];
        flint_mpn_zero(T.d, ext);
        T.n += ext; T.e -= ext; off = 0;
    }
    {
        /* materialize x << b into the tail of scratch */
        nn_ptr xs = scratch + x->n + (2 * (slong) i + 2) / FLINT_BITS + 6;
        slong xn = x->n;
        flint_mpn_copyi(xs, x->d, xn);
        if (b)
        {
            cy = mpn_lshift(xs, xs, xn, b);
            if (cy) xs[xn++] = cy;
        }
        while (off + xn > T.n)
            T.d[T.n++] = 0;
        cy = mpn_add(T.d + off, T.d + off, T.n - off, xs, xn);
        if (cy) T.d[T.n++] = cy;
    }

    /* T += x */
    off = x->e - T.e;
    if (off < 0)
    {
        slong ext = -off, t;
        for (t = T.n - 1; t >= 0; t--)
            T.d[t + ext] = T.d[t];
        flint_mpn_zero(T.d, ext);
        T.n += ext; T.e -= ext; off = 0;
    }
    while (off + x->n > T.n)
        T.d[T.n++] = 0;
    cy = mpn_add(T.d + off, T.d + off, T.n - off, x->d, x->n);
    if (cy) T.d[T.n++] = cy;

    flint_mpn_copyi(x->d, T.d, T.n);
    x->n = T.n; x->e = T.e;
    fx_norm(x);
    fx_trunc(x, cap, round_up);
}

static void
atanh_bs(fx_t * N, fx_t * DEN, fx_t * D, slong a, slong b, ulong i,
    slong cap, nn_ptr scratch, slong scratch_alloc)
{
    if (b - a <= TAB_BS_BLK)
    {
        slong j;
        slong bufn = cap + BUFSLACK(2 * (slong) i + 2);
        nn_ptr sc = scratch;

        /* W(j) = 1/d_j + q^-2 W(j+1):
           with V = DEN q^2:  N <- d_j N + V,  DEN <- d_j V.
           DEN stays EXACT within the block (it feeds both sides). */
        j = b - 1;
        N->d[0] = 1; N->n = 1; N->e = 0;
        DEN->d[0] = 2 * (ulong) j + 1; DEN->n = 1; DEN->e = 0;
        D->d[0] = 2 * (ulong) j + 1; D->n = 1; D->e = 0;

        for (j = b - 2; j >= a; j--)
        {
            ulong oj = 2 * (ulong) j + 1;
            ulong cy;
            slong off;

            /* DEN <- DEN q^2 (becomes V) */
            fx_mul_q2(DEN, i, bufn, 0, sc);

            /* N <- o_j N + V */
            cy = mpn_mul_1(N->d, N->d, N->n, oj);
            if (cy) N->d[N->n++] = cy;

            off = DEN->e - N->e;
            if (off < 0)
            {
                slong ext = -off, t;
                for (t = N->n - 1; t >= 0; t--)
                    N->d[t + ext] = N->d[t];
                flint_mpn_zero(N->d, ext);
                N->n += ext; N->e -= ext;
                off = 0;
            }
            while (off + DEN->n > N->n)
                N->d[N->n++] = 0;
            cy = mpn_add(N->d + off, N->d + off, N->n - off,
                DEN->d, DEN->n);
            if (cy) N->d[N->n++] = cy;
            fx_norm(N);
            fx_trunc(N, cap, 0);

            /* DEN <- o_j V (exact within the block) */
            cy = mpn_mul_1(DEN->d, DEN->d, DEN->n, oj);
            if (cy) DEN->d[DEN->n++] = cy;

            /* D <- o_j D */
            cy = mpn_mul_1(D->d, D->d, D->n, oj);
            if (cy) D->d[D->n++] = cy;
            fx_trunc(D, cap, 0);
        }

        /* one more q^2 converts to the per-term convention
           DEN = prod d_j q^(2(b-a)) */
        fx_mul_q2(DEN, i, bufn, 0, sc);
        fx_trunc(DEN, cap, 0);
        return;
    }

    {
        slong bufn = cap + BUFSLACK(2 * (slong) i + 2);
        slong m = a + (((b - a) / 2 + 1) & ~(slong) 1);
        fx_t N2, DEN2, D2, T1, T2;
        nn_ptr sc = scratch;

        N2.d = sc; sc += bufn;
        DEN2.d = sc; sc += bufn;
        D2.d = sc; sc += bufn;
        T1.d = sc; sc += 2 * bufn;
        T2.d = sc; sc += 2 * bufn;

        atanh_bs(N, DEN, D, a, m, i, cap, sc,
            scratch_alloc - (sc - scratch));
        atanh_bs(&N2, &DEN2, &D2, m, b, i, cap, sc,
            scratch_alloc - (sc - scratch));

        /* N = N1 DEN2 + N2 D1 */
        fx_mul(&T1, N, &DEN2, cap, 0);
        fx_mul(&T2, &N2, D, cap, 0);
        fx_add(N, &T1, &T2, cap);

        /* DEN = DEN1 DEN2 */
        fx_mul(&T1, DEN, &DEN2, cap, 0);
        flint_mpn_copyi(DEN->d, T1.d, T1.n);
        DEN->n = T1.n; DEN->e = T1.e;

        /* D = D1 D2 */
        fx_mul(&T1, D, &D2, cap, 0);
        flint_mpn_copyi(D->d, T1.d, T1.n);
        D->n = T1.n; D->e = T1.e;
    }
}

void
fixed_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
{
    slong prec = FLINT_BITS * n + 8;

    FLINT_ASSERT(i >= 1);
    FLINT_ASSERT(n >= 1);

    if (prec >= 6000 && i >= 30)
    {
        /* direct series: log(1 + 2^-i)
           = 2^-i sum_{a>=0} (-1)^a 2^(-i a) / (a + 1) */
        dyadic_t par;
        slong Nterms, F, depth, alloc, bufn;
        fx_t N, D;
        nn_ptr buf, scratch;

        Nterms = prec / (slong) i + 2;
        Nterms = (Nterms + 1) & ~(slong) 1;

        par.e = (slong) i;
        par.dbase = 1;
        par.dstep = 1;                        /* d_j = j + 1 */
        par.cap = n + 3;
        bufn = par.cap + BUFSLACK(par.e);

        for (depth = 1; (TAB_BS_BLK << depth) < Nterms; depth++)
            ;
        alloc = (depth + 2) * (6 * bufn + 40) + 4 * bufn + 40;

        buf = flint_malloc((2 * bufn + alloc) * sizeof(ulong));
        N.d = buf;
        D.d = buf + bufn;
        scratch = buf + 2 * bufn;

        dyadic_bs(&N, &D, &F, 0, Nterms, &par, scratch, alloc);
        dyadic_finish(res, n, &N, &D, F, (slong) i);

        flint_free(buf);
    }
    else
    {
        /* log(1 + 2^-i) = 2 atanh(1/q), q = 2^(i+1) + 1 */
        slong Nterms, depth, alloc, cap, bufn;
        fx_t N, DEN, D;
        nn_ptr buf, scratch;

        Nterms = prec / (2 * ((slong) i + 1)) + 2;
        Nterms = (Nterms + 1) & ~(slong) 1;

        cap = n + 3;
        bufn = cap + BUFSLACK(2 * (slong) i + 2);

        for (depth = 1; (TAB_BS_BLK << depth) < Nterms; depth++)
            ;
        alloc = (depth + 2) * (7 * bufn + 40) + 8 * bufn + 80;

        buf = flint_malloc((3 * bufn + alloc) * sizeof(ulong));
        N.d = buf;
        DEN.d = buf + bufn;
        D.d = buf + 2 * bufn;
        scratch = buf + 3 * bufn;

        atanh_bs(&N, &DEN, &D, 0, Nterms, i, cap, scratch, alloc);

        /* value = 2 atanh(1/q) = (2/q) sum q^(-2j)/d_j = 2 q W
           = 2 (N q) / DEN: the final q lands on the NUMERATOR
           (mpn_mul_1 when q fits a limb, else shift-and-add,
           rounding down).  DEN is then inflated past its
           rounded-down history -- at most a few events per node and
           per basecase step, each losing up to cap + 2 ulps now that
           the merges may go through the lower-approximating middle
           product -- so the denominator certainly reaches the true
           value and the quotient stays one-sided. */
        if (i + 1 < FLINT_BITS - 1)
        {
            ulong cy;
            cy = mpn_mul_1(N.d, N.d, N.n, (UWORD(1) << (i + 1)) + 1);
            if (cy)
                N.d[N.n++] = cy;
            fx_trunc(&N, cap + 2, 0);
        }
        else
        {
            fx_t T;
            slong off, t;
            ulong cy;

            T.d = scratch;
            flint_mpn_copyi(T.d, N.d, N.n);
            T.n = N.n; T.e = N.e;
            fx_shift_up(&T, (slong) i + 1);
            off = N.e - T.e;
            if (off < 0)
            {
                slong ext = -off;
                for (t = T.n - 1; t >= 0; t--)
                    T.d[t + ext] = T.d[t];
                flint_mpn_zero(T.d, ext);
                T.n += ext; T.e -= ext; off = 0;
            }
            while (off + N.n > T.n)
                T.d[T.n++] = 0;
            cy = mpn_add(T.d + off, T.d + off, T.n - off, N.d, N.n);
            if (cy) T.d[T.n++] = cy;
            flint_mpn_copyi(N.d, T.d, T.n);
            N.n = T.n; N.e = T.e;
            fx_norm(&N);
            fx_trunc(&N, cap + 2, 0);
        }

        {
            ulong cy;
            cy = mpn_add_1(DEN.d, DEN.d, DEN.n,
                (4 * (ulong) Nterms + 64) * ((ulong) cap + 4));
            if (cy)
                DEN.d[DEN.n++] = cy;
            fx_trunc(&DEN, cap + 2, 1);
        }

        dyadic_finish(res, n, &N, &DEN, 0, -1);

        flint_free(buf);
    }
}
