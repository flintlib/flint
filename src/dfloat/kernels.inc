/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The arithmetic kernels, written against an abstract element type so
   that they can be instantiated for scalars (KV = double) and for SIMD
   vectors (one lane per element).  Before including, define:

     KV             the element type
     KN(name)       name prefix for the functions
     K_ADD(a, b), K_SUB(a, b), K_MUL(a, b)
     K_FMS(a, b, c) a b - c, a single rounding
     K_NEG(a), K_ABS(a), K_ZERO, K_TWICE(a)
     K_SET_D(d)     the double d in every lane
     K_BOUND_BUMP_IF(a)  K_BOUND_BUMP where a is nonzero, else zero
     K_ANY_ZERO(a)  nonzero (int) if some lane of a is zero
     K_ANY_GAP(a, b)  nonzero (int) if some lane has a zero in a and a
                    nonzero in b (an internal zero of an expansion)

     KM             the number of components (the instantiation may
                    be for fewer components than the type's N)

   All kernels are branch-free except through K_ANY_ZERO and
   K_ANY_GAP, which are returned as a flag for the caller to handle
   (the removal of zero components before nonzero ones from the
   renormalized output, needed only after exact cancellation of an
   order sum: such a zero would let a later component be as large as
   u times the head, and the tapered kernels take the components at
   their positions).

   Every kernel computes its result with error-free transformations,
   and with track = 1 adds the magnitude of everything dropped to
   *err (exactly: nothing is added when nothing was dropped).  With
   track = 2 (bound mode, mul, sqr and mul_d only) the last order is
   accumulated with plain additions and their rounding errors are
   bounded by u times the accumulated values, which is cheaper but
   adds a nonzero bound whenever the last order has nonzero terms;
   see template.inc for the algorithms. */

#define K_TWO_SUM(s, e, a, b) \
    do { \
        KV _a = (a), _b = (b), _bb; \
        (s) = K_ADD(_a, _b); \
        _bb = K_SUB((s), _a); \
        (e) = K_ADD(K_SUB(_a, K_SUB((s), _bb)), K_SUB(_b, _bb)); \
    } while (0)

#define K_TWO_PROD(p, e, a, b) \
    do { \
        KV _a = (a), _b = (b); \
        (p) = K_MUL(_a, _b); \
        (e) = K_FMS(_a, _b, (p)); \
    } while (0)

#define K_ORDER_ACCUMULATE(acc, term, r, nr) \
    do { \
        KV _e; \
        K_TWO_SUM(acc, _e, acc, term); \
        (r)[(nr)++] = _e; \
    } while (0)

/* the accumulation of the last order in bound mode (track = 2): a
   plain addition, whose rounding error is at most u |acc| (for a
   normal acc; the bump for subnormal ones is added by the caller) */
#define K_ORDER_ACCUMULATE_BOUND(acc, term, bnd) \
    do { \
        (acc) = K_ADD(acc, term); \
        (bnd) = K_ADD(bnd, K_ABS(acc)); \
    } while (0)

#define K_ACCUM(k, acc, term, r, nr, bnd) \
    do { \
        if (track == 2 && (k) == KM - 1) \
            K_ORDER_ACCUMULATE_BOUND(acc, term, bnd); \
        else \
            K_ORDER_ACCUMULATE(acc, term, r, nr); \
    } while (0)

#define K_U 0x1p-53
#define K_BOUND_BUMP 0x1p-1069

/* Renormalize KM order sums S[0..KM-1] (|S[k]| a small multiple of
   u^k |S[0]|, so overlapping by at most a few bits) with a bottom-up
   VecSum: the result is weakly nonoverlapping (each t[i+1] is the
   rounding error of a sum whose rounded value is t[i], so
   |t[i+1]| <= ulp(t[i])/2, but a cancellation in t[i+1] can leave
   t[i+2] larger than t[i+1]), which is all the kernels need.  Nothing
   is dropped.  Returns nonzero if some head is zero, which happens
   only through cancellation between S[0] and the tail (the caller
   then shifts the components up).

   With strong set, a top-down VecSumErr pass follows, after which the
   expansion is ulp-nonoverlapping: |t[i+1]| <= ulp(t[i]) for every i,
   with zeros only at the end (TwoSum(0, b) = (b, 0) pushes zeros to
   the end).  TwoSum rather than Fast2Sum keeps the pass exact even
   in the near-cancellation cases where the Fast2Sum precondition
   could fail. */
FLINT_FORCE_INLINE int
KN(_renorm_orders)(KV * res, const KV * S, int strong)
{
    KV t[KM];
    int i;

    UNROLL
    for (i = 0; i < KM; i++)
        t[i] = S[i];
    UNROLL
    for (i = KM - 2; i >= 0; i--)
        K_TWO_SUM(t[i], t[i + 1], t[i], t[i + 1]);
    if (strong)
    {
        UNROLL
        for (i = 0; i < KM - 1; i++)
            K_TWO_SUM(t[i], t[i + 1], t[i], t[i + 1]);
    }
    UNROLL
    for (i = 0; i < KM; i++)
        res[i] = t[i];

    if (KM > 1)
    {
        int flag = K_ANY_ZERO(t[0]);
        UNROLL
        for (i = 1; i + 1 < KM; i++)
            flag |= K_ANY_GAP(t[i], t[i + 1]);
        return flag;
    }
    return 0;
}

/* res = x + y (or x - y with negate set): componentwise TwoSums give
   the sums s[k] and residuals e[k]; order k collects s[k], e[k-1] and
   the residuals of the order k-1 sum.  Dropped: the residuals of the
   last order and e[KM-1]. */
FLINT_FORCE_INLINE int
KN(_addsub_kernel)(KV * res, KV * err, const KV * x, const KV * y, int negate, int track, int strong)
{
    if (KM == 1)
    {
        KV s, e;
        K_TWO_SUM(s, e, x[0], negate ? K_NEG(y[0]) : y[0]);
        res[0] = s;
        if (track)
            *err = K_ADD(*err, K_ABS(e));
        return 0;
    }
    else
    {
        KV s[KM], e[KM], S[KM], acc, ee;
        KV R0[KM], R1[KM];
        KV * rprev = R0, * rnext = R1, * tmp;
        int i, k, nprev = 0, nnext;

        UNROLL
        for (i = 0; i < KM; i++)
            K_TWO_SUM(s[i], e[i], x[i], negate ? K_NEG(y[i]) : y[i]);

        S[0] = s[0];
        UNROLL
        for (k = 1; k < KM; k++)
        {
            nnext = 0;
            acc = s[k];
            K_ORDER_ACCUMULATE(acc, e[k - 1], rnext, nnext);
            UNROLL
            for (i = 0; i < nprev; i++)
                K_ORDER_ACCUMULATE(acc, rprev[i], rnext, nnext);
            S[k] = acc;
            tmp = rprev; rprev = rnext; rnext = tmp;
            nprev = nnext;
        }

        if (track)
        {
            ee = K_ABS(e[KM - 1]);
            UNROLL
            for (i = 0; i < nprev; i++)
                ee = K_ADD(ee, K_ABS(rprev[i]));
            *err = K_ADD(*err, ee);
        }

        return KN(_renorm_orders)(res, S, strong);
    }
}

/* res = x + y + z: like the sum, with two residuals e1[k], e2[k] per
   component from the two TwoSums; order k collects s[k], e1[k-1],
   e2[k-1] and the residuals of order k-1.  Dropped: the residuals of
   the last order and e1[KM-1], e2[KM-1]. */
FLINT_FORCE_INLINE int
KN(_add3_kernel)(KV * res, KV * err, const KV * x, const KV * y, const KV * z, int track, int strong)
{
    if (KM == 1)
    {
        KV s, e1, e2;
        K_TWO_SUM(s, e1, x[0], y[0]);
        K_TWO_SUM(s, e2, s, z[0]);
        res[0] = s;
        if (track)
            *err = K_ADD(*err, K_ADD(K_ABS(e1), K_ABS(e2)));
        return 0;
    }
    else
    {
        KV s[KM], e1[KM], e2[KM], S[KM], acc, ee;
        KV R0[2 * KM], R1[2 * KM];
        KV * rprev = R0, * rnext = R1, * tmp;
        int i, k, nprev = 0, nnext;

        UNROLL
        for (i = 0; i < KM; i++)
        {
            K_TWO_SUM(s[i], e1[i], x[i], y[i]);
            K_TWO_SUM(s[i], e2[i], s[i], z[i]);
        }

        S[0] = s[0];
        UNROLL
        for (k = 1; k < KM; k++)
        {
            nnext = 0;
            acc = s[k];
            K_ORDER_ACCUMULATE(acc, e1[k - 1], rnext, nnext);
            K_ORDER_ACCUMULATE(acc, e2[k - 1], rnext, nnext);
            UNROLL
            for (i = 0; i < nprev; i++)
                K_ORDER_ACCUMULATE(acc, rprev[i], rnext, nnext);
            S[k] = acc;
            tmp = rprev; rprev = rnext; rnext = tmp;
            nprev = nnext;
        }

        if (track)
        {
            ee = K_ADD(K_ABS(e1[KM - 1]), K_ABS(e2[KM - 1]));
            UNROLL
            for (i = 0; i < nprev; i++)
                ee = K_ADD(ee, K_ABS(rprev[i]));
            *err = K_ADD(*err, ee);
        }

        return KN(_renorm_orders)(res, S, strong);
    }
}

/* res = xs[0] + ... + xs[nx-1] (nx >= 2): componentwise, nx - 1
   TwoSums per component whose residuals all go to the next order.
   Dropped: the residuals of the last order and the residuals of the
   component KM-1. */
FLINT_FORCE_INLINE int
KN(_addn_kernel)(KV * res, KV * err, const KV * const * xs, int nx, int track, int strong)
{
    KV s[KM], ec[KM][8], S[KM], acc, ee;
    KV R0[8 * KM], R1[8 * KM];
    KV * rprev = R0, * rnext = R1, * tmp;
    int i, k, m, nprev = 0, nnext;

    UNROLL
    for (i = 0; i < KM; i++)
    {
        s[i] = xs[0][i];
        UNROLL
        for (m = 1; m < nx; m++)
            K_TWO_SUM(s[i], ec[i][m - 1], s[i], xs[m][i]);
    }

    if (KM == 1)
    {
        res[0] = s[0];
        if (track)
        {
            ee = K_ZERO;
            for (m = 1; m < nx; m++)
                ee = K_ADD(ee, K_ABS(ec[0][m - 1]));
            *err = K_ADD(*err, ee);
        }
        return 0;
    }

    S[0] = s[0];
    UNROLL
    for (k = 1; k < KM; k++)
    {
        nnext = 0;
        acc = s[k];
        UNROLL
        for (m = 1; m < nx; m++)
            K_ORDER_ACCUMULATE(acc, ec[k - 1][m - 1], rnext, nnext);
        UNROLL
        for (i = 0; i < nprev; i++)
            K_ORDER_ACCUMULATE(acc, rprev[i], rnext, nnext);
        S[k] = acc;
        tmp = rprev; rprev = rnext; rnext = tmp;
        nprev = nnext;
    }

    if (track)
    {
        ee = K_ZERO;
        UNROLL
        for (m = 1; m < nx; m++)
            ee = K_ADD(ee, K_ABS(ec[KM - 1][m - 1]));
        UNROLL
        for (i = 0; i < nprev; i++)
            ee = K_ADD(ee, K_ABS(rprev[i]));
        *err = K_ADD(*err, ee);
    }

    return KN(_renorm_orders)(res, S, strong);
}

/* res = x y.  The partial products x[i] y[j] of order i + j < KM are
   computed exactly (TwoProd) and summed order by order with TwoSums;
   the residuals of order k feed the sum of order k + 1, so that only
   the residuals of the last order, the low parts of the order KM - 1
   products, and the products of order >= KM (bounded by their rounded
   magnitudes) are dropped.  Requires that no partial product
   underflows (all nonzero components >= 2^-480 in absolute value)
   when track is set. */
FLINT_FORCE_INLINE int
KN(_mul_kernel)(KV * res, KV * err, const KV * x, const KV * y, int track, int strong)
{
    if (KM == 1)
    {
        KV p, e;
        K_TWO_PROD(p, e, x[0], y[0]);
        res[0] = p;
        if (track)
            *err = K_ADD(*err, K_ABS(e));
        return 0;
    }
    else
    {
        KV h[KM][KM], l[KM][KM];
        KV S[KM], acc, e, bnd;
        KV R0[KM * KM], R1[KM * KM];
        KV * rprev = R0, * rnext = R1, * tmp;
        int i, j, k, nprev = 0, nnext;

        UNROLL
        for (i = 0; i < KM; i++)
        {
            UNROLL
            for (j = 0; j + i < KM; j++)
                K_TWO_PROD(h[i][j], l[i][j], x[i], y[j]);
        }

        S[0] = h[0][0];

        bnd = K_ZERO;
        UNROLL
        for (k = 1; k < KM; k++)
        {
            nnext = 0;
            acc = h[0][k];
            UNROLL
            for (i = 1; i <= k; i++)
                K_ACCUM(k, acc, h[i][k - i], rnext, nnext, bnd);
            UNROLL
            for (i = 0; i < k; i++)
                K_ACCUM(k, acc, l[i][k - 1 - i], rnext, nnext, bnd);
            UNROLL
            for (i = 0; i < nprev; i++)
                K_ACCUM(k, acc, rprev[i], rnext, nnext, bnd);
            S[k] = acc;
            tmp = rprev; rprev = rnext; rnext = tmp;
            nprev = nnext;
        }

        if (track)
        {
            if (track == 2)
                e = K_ADD(K_MUL(bnd, K_SET_D(K_U)), K_BOUND_BUMP_IF(bnd));
            else
            {
                e = K_ZERO;
                UNROLL
                for (i = 0; i < nprev; i++)
                    e = K_ADD(e, K_ABS(rprev[i]));
            }
            UNROLL
            for (i = 0; i < KM; i++)
                e = K_ADD(e, K_ABS(l[i][KM - 1 - i]));
            UNROLL
            for (i = 1; i < KM; i++)
            {
                UNROLL
                for (j = KM - i; j < KM; j++)
                    e = K_ADD(e, K_ABS(K_MUL(x[i], y[j])));
            }
            *err = K_ADD(*err, e);
        }

        return KN(_renorm_orders)(res, S, strong);
    }
}

/* res = x^2: the same with the symmetric products doubled (exact). */
FLINT_FORCE_INLINE int
KN(_sqr_kernel)(KV * res, KV * err, const KV * x, int track, int strong)
{
    if (KM == 1)
    {
        KV p, e;
        K_TWO_PROD(p, e, x[0], x[0]);
        res[0] = p;
        if (track)
            *err = K_ADD(*err, K_ABS(e));
        return 0;
    }
    else
    {
        KV h[KM][KM], l[KM][KM];
        KV S[KM], acc, e, bnd;
        KV R0[KM * KM], R1[KM * KM];
        KV * rprev = R0, * rnext = R1, * tmp;
        int i, j, k, nprev = 0, nnext;

        UNROLL
        for (i = 0; i < KM; i++)
        {
            UNROLL
            for (j = i; j + i < KM; j++)
            {
                K_TWO_PROD(h[i][j], l[i][j], x[i], x[j]);
                if (i != j)
                {
                    h[i][j] = K_TWICE(h[i][j]);
                    l[i][j] = K_TWICE(l[i][j]);
                }
            }
        }

        S[0] = h[0][0];

        bnd = K_ZERO;
        UNROLL
        for (k = 1; k < KM; k++)
        {
            nnext = 0;
            acc = h[0][k];
            UNROLL
            for (i = 1; 2 * i <= k; i++)
                K_ACCUM(k, acc, h[i][k - i], rnext, nnext, bnd);
            UNROLL
            for (i = 0; 2 * i <= k - 1; i++)
                K_ACCUM(k, acc, l[i][k - 1 - i], rnext, nnext, bnd);
            UNROLL
            for (i = 0; i < nprev; i++)
                K_ACCUM(k, acc, rprev[i], rnext, nnext, bnd);
            S[k] = acc;
            tmp = rprev; rprev = rnext; rnext = tmp;
            nprev = nnext;
        }

        if (track)
        {
            if (track == 2)
                e = K_ADD(K_MUL(bnd, K_SET_D(K_U)), K_BOUND_BUMP_IF(bnd));
            else
            {
                e = K_ZERO;
                UNROLL
                for (i = 0; i < nprev; i++)
                    e = K_ADD(e, K_ABS(rprev[i]));
            }
            UNROLL
            for (i = 0; 2 * i <= KM - 1; i++)
                e = K_ADD(e, K_ABS(l[i][KM - 1 - i]));
            UNROLL
            for (i = 1; i < KM; i++)
            {
                UNROLL
                for (j = FLINT_MAX(i, KM - i); j < KM; j++)
                {
                    KV _p = K_ABS(K_MUL(x[i], x[j]));
                    e = K_ADD(e, (i == j) ? _p : K_TWICE(_p));
                }
            }
            *err = K_ADD(*err, e);
        }

        return KN(_renorm_orders)(res, S, strong);
    }
}

/* res = x y for a scalar y: exact products h[k] + l[k]; order k
   collects h[k], l[k-1] and the previous residuals.  Dropped: the last
   residuals and l[KM-1].  Requires no underflow when track is set. */
FLINT_FORCE_INLINE int
KN(_mul_d_kernel)(KV * res, KV * err, const KV * x, KV y, int track, int strong)
{
    if (KM == 1)
    {
        KV p, e;
        K_TWO_PROD(p, e, x[0], y);
        res[0] = p;
        if (track)
            *err = K_ADD(*err, K_ABS(e));
        return 0;
    }
    else
    {
        KV h[KM], l[KM], S[KM], acc, ee;
        KV R0[KM], R1[KM];
        KV * rprev = R0, * rnext = R1, * tmp;
        int i, k, nprev = 0, nnext;

        UNROLL
        for (i = 0; i < KM; i++)
            K_TWO_PROD(h[i], l[i], x[i], y);

        S[0] = h[0];
        UNROLL
        for (k = 1; k < KM; k++)
        {
            nnext = 0;
            acc = h[k];
            K_ORDER_ACCUMULATE(acc, l[k - 1], rnext, nnext);
            UNROLL
            for (i = 0; i < nprev; i++)
                K_ORDER_ACCUMULATE(acc, rprev[i], rnext, nnext);
            S[k] = acc;
            tmp = rprev; rprev = rnext; rnext = tmp;
            nprev = nnext;
        }

        if (track)
        {
            ee = K_ABS(l[KM - 1]);
            UNROLL
            for (i = 0; i < nprev; i++)
                ee = K_ADD(ee, K_ABS(rprev[i]));
            *err = K_ADD(*err, ee);
        }

        return KN(_renorm_orders)(res, S, strong);
    }
}

#undef K_TWO_SUM
#undef K_TWO_PROD
#undef K_ORDER_ACCUMULATE
#undef K_ORDER_ACCUMULATE_BOUND
#undef K_ACCUM
#undef K_U
#undef K_BOUND_BUMP
