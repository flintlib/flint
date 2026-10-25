/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "padic_radix.h"
#include "impl.h"

/*
    Min-plus convolutions of exponent sequences, for the error propagation
    of polynomial products: c_k = sum_i a_i b_(k-i) has error term
    min_i (N(a_i) + v(b_(k-i))) (and the symmetric term), a min-plus
    convolution. Entries >= PADIC_RADIX_PACK_NONE stand for +infinity
    (no term), and results >= PACK_NONE mean that there is no finite term.

    The exact convolution costs O(alen blen). The lower bound replaces each
    sequence by its lower convex hull (the largest convex minorant, defined
    on the interval spanned by the finite entries): the min-plus convolution
    of two convex functions is again convex, with the edges of both hulls
    merged by slope (a Minkowski sum of the epigraphs), so it costs
    O(alen + blen), and it is a lower bound since the hulls are minorants.
    It is exact when both sequences are convex (in particular constant or
    linear, e.g. uniform valuations and precisions). Rounded up to integers,
    it remains a valid lower bound for the integer exponents.
*/

/* sign of (y1 * x2 - y2 * x1), exactly */
static int
_cross_sign(slong y1, slong x1, slong y2, slong x2)
{
    ulong h1, l1, h2, l2;
    smul_ppmm(h1, l1, y1, x2);
    smul_ppmm(h2, l2, y2, x1);
    if (h1 != h2)
        return ((slong) h1 < (slong) h2) ? -1 : 1;
    if (l1 != l2)
        return (l1 < l2) ? -1 : 1;
    return 0;
}

/* lower convex hull of the finite points (i, y[i]); vertices in hx, hy
   (increasing x); returns the number of vertices */
static slong
_lower_hull(slong * hx, slong * hy, const slong * y, slong len)
{
    slong i, m = 0;

    for (i = 0; i < len; i++)
    {
        if (y[i] >= PADIC_RADIX_PACK_NONE)
            continue;

        /* pop while the last two vertices and (i, y[i]) do not turn left,
           i.e. slope(h[m-2], h[m-1]) >= slope(h[m-1], i) */
        while (m >= 2 && _cross_sign(hy[m - 1] - hy[m - 2], hx[m - 1] - hx[m - 2],
                                     y[i] - hy[m - 1], i - hx[m - 1]) >= 0)
            m--;

        hx[m] = i;
        hy[m] = y[i];
        m++;
    }

    return m;
}

/* ceil(dy t / dx) for dx > 0, 0 <= t <= dx, without overflow (|dy| below
   2^62 and dx below 2^31): with dy = q dx + r, 0 <= r < dx, it is
   q t + ceil(r t / dx), where r t < dx^2 */
static slong
_interp_ceil(slong dy, slong t, slong dx)
{
    slong q = dy / dx, r = dy - q * dx;
    if (r < 0)
    {
        q--;
        r += dx;
    }
    return q * t + (r * t + dx - 1) / dx;
}

void
_padic_radix_minplus_conv_lower(slong * res, const slong * a, slong alen,
    const slong * b, slong blen, slong klo, slong khi)
{
    slong * ax, * ay, * bx, * by;
    slong am, bm, k, ia, ib, x, y;
    TMP_INIT;

    TMP_START;
    ax = TMP_ALLOC(sizeof(slong) * 2 * (alen + blen));
    ay = ax + alen;
    bx = ay + alen;
    by = bx + blen;

    am = _lower_hull(ax, ay, a, alen);
    bm = _lower_hull(bx, by, b, blen);

    if (am == 0 || bm == 0)
    {
        for (k = klo; k < khi; k++)
            res[k - klo] = 2 * PADIC_RADIX_PACK_BIG;
        TMP_END;
        return;
    }

    /* walk the merged edges from (ax[0] + bx[0], ay[0] + by[0]) */
    x = ax[0] + bx[0];
    y = ay[0] + by[0];
    ia = 0;
    ib = 0;
    k = klo;

    for (; k < khi && k < x; k++)
        res[k - klo] = 2 * PADIC_RADIX_PACK_BIG;

    while (k < khi)
    {
        slong dx, dy;

        if (k == x)
            res[k - klo] = y;

        if (ia == am - 1 && ib == bm - 1)
            break;

        /* next edge: the one of smaller slope */
        if (ib == bm - 1 || (ia < am - 1 &&
            _cross_sign(ay[ia + 1] - ay[ia], ax[ia + 1] - ax[ia],
                        by[ib + 1] - by[ib], bx[ib + 1] - bx[ib]) <= 0))
        {
            dx = ax[ia + 1] - ax[ia];
            dy = ay[ia + 1] - ay[ia];
            ia++;
        }
        else
        {
            dx = bx[ib + 1] - bx[ib];
            dy = by[ib + 1] - by[ib];
            ib++;
        }

        /* integer points strictly inside the edge, then its end */
        for (k = FLINT_MAX(k, x + 1); k < khi && k < x + dx; k++)
            res[k - klo] = y + _interp_ceil(dy, k - x, dx);

        x += dx;
        y += dy;
    }

    for (k = FLINT_MAX(k, x + 1); k < khi; k++)
        res[k - klo] = 2 * PADIC_RADIX_PACK_BIG;

    TMP_END;
}

void
_padic_radix_minplus_conv_exact(slong * res, const slong * a, slong alen,
    const slong * b, slong blen, slong klo, slong khi)
{
    slong k, i;

    for (k = klo; k < khi; k++)
    {
        slong ilo = FLINT_MAX(0, k - blen + 1), ihi = FLINT_MIN(k, alen - 1);
        slong m = 2 * PADIC_RADIX_PACK_BIG;

        for (i = ilo; i <= ihi; i++)
            m = FLINT_MIN(m, a[i] + b[k - i]);

        res[k - klo] = m;
    }
}

/*
    Exact min-plus convolution of sequences whose finite values span a small
    range R (after subtracting the minima): with F = sum_i x^i y^(a_i) and
    G = sum_j x^j y^(b_j), the minimum over i of a_i + b_(k-i) is the lowest
    power of y in the coefficient of x^k of F G. Bivariate Kronecker
    substitution with fields of fbits bits for the powers of y (a count of
    pairs fits) and 2R + 1 fields per power of x turns this into one integer
    multiplication of O((alen + blen) R log(len)) bits.
*/
static void
_minplus_conv_ks(slong * res, const slong * a, slong alen, slong amin,
    const slong * b, slong blen, slong bmin, slong R, slong klo, slong khi)
{
    slong S = 2 * R + 1, fbits, an, bn, rn, i, k, t;
    nn_ptr A, Bp, C;
    TMP_INIT;

    fbits = FLINT_BIT_COUNT(FLINT_MIN(alen, blen)) + 1;
    an = (alen * S * fbits + FLINT_BITS - 1) / FLINT_BITS + 1;
    bn = (blen * S * fbits + FLINT_BITS - 1) / FLINT_BITS + 1;
    rn = an + bn;

    TMP_START;
    A = TMP_ALLOC(sizeof(ulong) * (an + bn + rn));
    Bp = A + an;
    C = Bp + bn;
    flint_mpn_zero(A, an + bn);

    for (i = 0; i < alen; i++)
        if (a[i] < PADIC_RADIX_PACK_NONE)
        {
            ulong bit = (ulong) (i * S + (a[i] - amin)) * fbits;
            A[bit / FLINT_BITS] |= UWORD(1) << (bit % FLINT_BITS);
        }

    for (i = 0; i < blen; i++)
        if (b[i] < PADIC_RADIX_PACK_NONE)
        {
            ulong bit = (ulong) (i * S + (b[i] - bmin)) * fbits;
            Bp[bit / FLINT_BITS] |= UWORD(1) << (bit % FLINT_BITS);
        }

    if (an >= bn)
        flint_mpn_mul(C, A, an, Bp, bn);
    else
        flint_mpn_mul(C, Bp, bn, A, an);

    for (k = klo; k < khi; k++)
    {
        slong m = 2 * PADIC_RADIX_PACK_BIG;

        for (t = 0; t < S; t++)
        {
            ulong bit = (ulong) (k * S + t) * fbits;
            ulong w = C[bit / FLINT_BITS] >> (bit % FLINT_BITS);
            if ((bit % FLINT_BITS) + fbits > FLINT_BITS)
                w |= C[bit / FLINT_BITS + 1] << (FLINT_BITS - (bit % FLINT_BITS));
            if (w & ((UWORD(1) << fbits) - 1))
            {
                m = amin + bmin + t;
                break;
            }
        }

        res[k - klo] = m;
    }

    TMP_END;
}

/* range of the finite entries of y[i] - lambda i; returns 0 if none */
static int
_detrended_range(slong * ymin, slong * ymax, const slong * y, slong len, slong lambda)
{
    slong i, lo = WORD_MAX, hi = WORD_MIN;
    for (i = 0; i < len; i++)
        if (y[i] < PADIC_RADIX_PACK_NONE)
        {
            slong z = y[i] - lambda * i;
            lo = FLINT_MIN(lo, z);
            hi = FLINT_MAX(hi, z);
        }
    *ymin = lo;
    *ymax = hi;
    return lo != WORD_MAX;
}

/* average slope of the finite entries (first to last), rounded down */
static slong
_trend(const slong * y, slong len)
{
    slong i0 = -1, i1 = -1, i;
    for (i = 0; i < len; i++)
        if (y[i] < PADIC_RADIX_PACK_NONE)
        {
            if (i0 < 0)
                i0 = i;
            i1 = i;
        }
    if (i0 < 0 || i1 == i0)
        return 0;
    {
        slong d = y[i1] - y[i0], q = d / (i1 - i0);
        if (q * (i1 - i0) > d)
            q--;
        return q;
    }
}

#define MINPLUS_SMALL_PAIRS 4096
#define MINPLUS_KS_MAX_RANGE 64
#define MINPLUS_EDGE_TERMS 32
#define MINPLUS_BANDS 32
#define MINPLUS_BAND_PAIRS 65536

#define INF (2 * PADIC_RADIX_PACK_BIG)

/* indices of the finite entries */
static slong
_finite_indices(slong * idx, const slong * y, slong len)
{
    slong i, n = 0;
    for (i = 0; i < len; i++)
        if (y[i] < PADIC_RADIX_PACK_NONE)
            idx[n++] = i;
    return n;
}

/* first position in the sorted list idx[0, n) with idx[pos] >= x */
static slong
_lower_index(const slong * idx, slong n, slong x)
{
    slong lo = 0, hi = n;
    while (lo < hi)
    {
        slong mid = lo + (hi - lo) / 2;
        if (idx[mid] < x)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

/* res[k - klo] = min(res[k - klo], a_i + b_j) over the listed finite
   entries with klo <= i + j < khi: O(na nb) */
static void
_minplus_sparse_update(slong * res, const slong * a, const slong * ai, slong na,
    const slong * b, const slong * bj, slong nb, slong klo, slong khi)
{
    slong s, t;

    for (s = 0; s < na; s++)
    {
        slong i = ai[s], x = a[i];
        slong t0 = _lower_index(bj, nb, klo - i);
        slong t1 = _lower_index(bj, nb, khi - i);

        for (t = t0; t < t1; t++)
        {
            slong k = i + bj[t] - klo, y = x + b[bj[t]];
            res[k] = FLINT_MIN(res[k], y);
        }
    }
}

/* whether the finite entries are contiguous and convex */
static int
_contiguous_convex(const slong * y, const slong * idx, slong n)
{
    slong t;

    if (n == 0 || idx[n - 1] - idx[0] != n - 1)
        return 0;

    for (t = 1; t + 1 < n; t++)
    {
        /* y[i-1] + y[i+1] >= 2 y[i] (values are moderate: no overflow) */
        if (y[idx[t - 1]] + y[idx[t + 1]] < 2 * y[idx[t]])
            return 0;
    }

    return 1;
}

/* rough costs in nanoseconds: the Kronecker substitution costs about as
   much per bit of its operands as the quadratic algorithm per pair */
static double
_ks_cost(slong alen, slong blen, slong nout, slong R)
{
    double fbits = FLINT_BIT_COUNT(FLINT_MIN(alen, blen)) + 1;
    double bits = (double) (alen + blen) * (2.0 * R + 1) * fbits;
    return 300.0 + 0.7 * bits + 0.1 * (double) nout * (2.0 * R + 1);
}

static double
_direct_cost(slong na, slong nb)
{
    return 50.0 + 0.7 * (double) na * (double) nb;
}

/*
    Exact min-plus convolution by bands of values: with the (detrended)
    finite values split into bands [lo + s W, lo + (s + 1) W), the result is
    the minimum over pairs of bands (s, t) of the convolutions restricted to
    those bands. Pairs are processed by increasing s + t; an output whose
    current value is at most lo_a + lo_b + (s + t) W cannot be improved by
    the remaining pairs and is settled, and a pair whose output range only
    contains settled outputs (or outputs without any finite term) is
    skipped. Each pair uses the cheaper of the direct sparse product and the
    Kronecker substitution (range W). For few distinct bands or clustered
    values this is close to one Kronecker substitution; in the worst case
    (values spread at random over a large range) it does about the work of
    the quadratic algorithm.
*/
static void
_minplus_conv_bands(slong * res, const slong * a, slong alen, slong amin, slong amax,
    const slong * b, slong blen, slong bmin, slong bmax, slong klo, slong khi)
{
    slong R = FLINT_MAX(amax - amin, bmax - bmin);
    slong W, Ba, Bb, nout = khi - klo, i, k, g, s;
    slong * acnt, * bcnt, * aoff, * boff, * alist, * blist, * afill, * bfill;
    slong * reach, * unset, * tmp, * sub, * alo, * ahi, * blo, * bhi;
    slong na, nb;
    TMP_INIT;

    W = (R + MINPLUS_BANDS) / MINPLUS_BANDS;
    Ba = (amax - amin) / W + 1;
    Bb = (bmax - bmin) / W + 1;

    TMP_START;

    acnt = TMP_ALLOC(sizeof(slong) * (2 * (Ba + 1) + 2 * (Bb + 1) + Ba + Bb));
    aoff = acnt + (Ba + 1);
    bcnt = aoff + (Ba + 1);
    boff = bcnt + (Bb + 1);
    afill = boff + (Bb + 1);
    bfill = afill + Ba;

    alo = TMP_ALLOC(sizeof(slong) * 2 * (Ba + Bb));
    ahi = alo + Ba;
    blo = ahi + Ba;
    bhi = blo + Bb;

    alist = TMP_ALLOC(sizeof(slong) * (alen + blen + 3 * nout + 2 + alen + blen));
    blist = alist + alen;
    reach = blist + blen;
    unset = reach + nout;        /* nout + 1 prefix sums */
    tmp = unset + nout + 1;
    sub = tmp + nout;            /* alen + blen */

    /* bucket the indices by band (counting sort keeps them increasing) */
    for (s = 0; s <= Ba; s++) acnt[s] = 0;
    for (s = 0; s <= Bb; s++) bcnt[s] = 0;
    for (i = 0; i < alen; i++)
        if (a[i] < PADIC_RADIX_PACK_NONE)
            acnt[(a[i] - amin) / W]++;
    for (i = 0; i < blen; i++)
        if (b[i] < PADIC_RADIX_PACK_NONE)
            bcnt[(b[i] - bmin) / W]++;
    aoff[0] = 0;
    for (s = 0; s < Ba; s++) aoff[s + 1] = aoff[s] + acnt[s];
    boff[0] = 0;
    for (s = 0; s < Bb; s++) boff[s + 1] = boff[s] + bcnt[s];
    na = aoff[Ba];
    nb = boff[Bb];
    for (s = 0; s < Ba; s++) afill[s] = aoff[s];
    for (s = 0; s < Bb; s++) bfill[s] = boff[s];
    for (s = 0; s < Ba; s++) { alo[s] = WORD_MAX; ahi[s] = WORD_MIN; }
    for (s = 0; s < Bb; s++) { blo[s] = WORD_MAX; bhi[s] = WORD_MIN; }
    for (i = 0; i < alen; i++)
        if (a[i] < PADIC_RADIX_PACK_NONE)
        {
            slong q = (a[i] - amin) / W;
            alist[afill[q]++] = i;
            alo[q] = FLINT_MIN(alo[q], a[i]);
            ahi[q] = FLINT_MAX(ahi[q], a[i]);
        }
    for (i = 0; i < blen; i++)
        if (b[i] < PADIC_RADIX_PACK_NONE)
        {
            slong q = (b[i] - bmin) / W;
            blist[bfill[q]++] = i;
            blo[q] = FLINT_MIN(blo[q], b[i]);
            bhi[q] = FLINT_MAX(bhi[q], b[i]);
        }

    /* outputs with at least one finite term: support convolution */
    for (i = 0; i < alen; i++)
        sub[i] = (a[i] < PADIC_RADIX_PACK_NONE) ? 0 : INF;
    for (i = 0; i < blen; i++)
        sub[alen + i] = (b[i] < PADIC_RADIX_PACK_NONE) ? 0 : INF;
    if (_direct_cost(na, nb) < _ks_cost(alen, blen, nout, 0))
    {
        slong * ai = TMP_ALLOC(sizeof(slong) * (na + nb)), * bj = ai + na;
        _finite_indices(ai, a, alen);
        _finite_indices(bj, b, blen);
        for (k = 0; k < nout; k++)
            reach[k] = INF;
        _minplus_sparse_update(reach, sub, ai, na, sub + alen, bj, nb, klo, khi);
    }
    else
        _minplus_conv_ks(reach, sub, alen, 0, sub + alen, blen, 0, 0, klo, khi);

    for (k = 0; k < nout; k++)
        res[k] = INF;

    for (g = 0; g <= Ba + Bb - 2; g++)
    {
        slong T = amin + bmin + g * W, total;

        /* unsettled outputs: reachable, current value above T */
        unset[0] = 0;
        for (k = 0; k < nout; k++)
            unset[k + 1] = unset[k] + (reach[k] < PADIC_RADIX_PACK_NONE && res[k] > T);
        total = unset[nout];
        if (total == 0)
            break;

        for (s = FLINT_MAX(0, g - (Bb - 1)); s <= FLINT_MIN(g, Ba - 1); s++)
        {
            slong t = g - s;
            slong ns = acnt[s], nt = bcnt[t];
            slong ilo, ihi, jlo, jhi, olo, ohi, Rst;

            if (ns == 0 || nt == 0)
                continue;

            ilo = alist[aoff[s]];
            ihi = alist[aoff[s] + ns - 1];
            jlo = blist[boff[t]];
            jhi = blist[boff[t] + nt - 1];

            olo = FLINT_MAX(klo, ilo + jlo);
            ohi = FLINT_MIN(khi, ihi + jhi + 1);
            if (olo >= ohi || unset[ohi - klo] == unset[olo - klo])
                continue;

            Rst = FLINT_MAX(ahi[s] - alo[s], bhi[t] - blo[t]);

            if (_direct_cost(ns, nt) <= _ks_cost(ihi - ilo + 1, jhi - jlo + 1, ohi - olo, Rst))
            {
                _minplus_sparse_update(res, a, alist + aoff[s], ns, b, blist + boff[t], nt, klo, khi);
            }
            else
            {
                slong la = ihi - ilo + 1, lb = jhi - jlo + 1, q;

                for (i = 0; i < la; i++)
                    sub[i] = INF;
                for (i = 0; i < lb; i++)
                    sub[la + i] = INF;
                for (q = aoff[s]; q < aoff[s] + ns; q++)
                    sub[alist[q] - ilo] = a[alist[q]];
                for (q = boff[t]; q < boff[t] + nt; q++)
                    sub[la + blist[q] - jlo] = b[blist[q]];

                _minplus_conv_ks(tmp, sub, la, alo[s], sub + la, lb, blo[t], Rst,
                    olo - ilo - jlo, ohi - ilo - jlo);

                for (k = olo; k < ohi; k++)
                    res[k - klo] = FLINT_MIN(res[k - klo], tmp[k - olo]);
            }
        }
    }

    TMP_END;
}

/* the common slope (among 0 and the average slopes of a and b) that
   minimizes the larger of the detrended ranges */
static slong
_best_trend(slong * amin, slong * amax, slong * bmin, slong * bmax,
    const slong * a, slong alen, const slong * b, slong blen)
{
    slong cand[3], best = 0, bestR = WORD_MAX, c;

    cand[0] = 0;
    cand[1] = _trend(a, alen);
    cand[2] = _trend(b, blen);

    for (c = 0; c < 3; c++)
    {
        slong l = cand[c], x0, x1, y0, y1, R;

        if (c > 0 && (l == 0 || FLINT_ABS(l) >= PADIC_RADIX_PACK_MODERATE / FLINT_MAX(alen + blen, 1)))
            continue;

        _detrended_range(&x0, &x1, a, alen, l);
        _detrended_range(&y0, &y1, b, blen, l);
        R = FLINT_MAX(x1 - x0, y1 - y0);

        if (c == 0 || R < bestR)
        {
            bestR = R;
            best = l;
            *amin = x0; *amax = x1; *bmin = y0; *bmax = y1;
        }
    }

    return best;
}

/*
    res[k - klo] = min_i (a_i + b_(k-i)) for klo <= k < khi.

    With exact != 0 the result is always exact: by the direct product of
    the finite entries when there are few pairs, by the convex hulls when
    both sequences are convex without gaps, by one Kronecker substitution
    when the values (after subtracting a common linear trend) span a small
    range, and otherwise by bands of values (above). The cost is
    O(n log n) in the regular cases and degrades towards the quadratic cost
    for irregular ones.

    With exact = 0 the result may be a lower bound, computed in quasi-linear
    time in the irregular cases: the larger of the convex hull bound and the
    Kronecker substitution on values rounded down to a multiple of delta
    (within 2 delta), with the outputs near the ends computed exactly; the
    outputs without any finite term are always recognised. Returns 1 if the
    result is known to be exact.
*/
int
_padic_radix_minplus_conv(slong * res, const slong * a, slong alen,
    const slong * b, slong blen, slong klo, slong khi, int exact)
{
    slong lambda, amin = 0, amax = 0, bmin = 0, bmax = 0, R, k, i, i0, j0, nout, na, nb;
    slong * ad, * bd, * ai, * bj;
    int result = 1;
    TMP_INIT;

    nout = khi - klo;
    if (nout <= 0)
        return 1;

    /* only a[i0, ...) and b[j0, ...) reach the outputs */
    i0 = FLINT_MAX(0, klo - (blen - 1));
    j0 = FLINT_MAX(0, klo - (alen - 1));
    alen = FLINT_MIN(alen, khi) - i0;
    blen = FLINT_MIN(blen, khi) - j0;

    if (alen <= 0 || blen <= 0)
    {
        for (k = 0; k < nout; k++)
            res[k] = INF;
        return 1;
    }

    a += i0;
    b += j0;
    klo -= i0 + j0;
    khi -= i0 + j0;

    TMP_START;

    ai = TMP_ALLOC(sizeof(slong) * (2 * (alen + blen) + nout));
    bj = ai + alen;
    ad = bj + blen;
    bd = ad + alen;

    na = _finite_indices(ai, a, alen);
    nb = _finite_indices(bj, b, blen);

    if (na == 0 || nb == 0)
    {
        for (k = 0; k < nout; k++)
            res[k] = INF;
        goto cleanup;
    }

    /* few pairs */
    if ((double) na * (double) nb <= MINPLUS_SMALL_PAIRS || nout <= 2)
    {
        if (na == alen && nb == blen)
            _padic_radix_minplus_conv_exact(res, a, alen, b, blen, klo, khi);
        else
        {
            for (k = 0; k < nout; k++)
                res[k] = INF;
            _minplus_sparse_update(res, a, ai, na, b, bj, nb, klo, khi);
        }
        goto cleanup;
    }

    /* convex sequences without gaps: the hull walk is exact */
    if (_contiguous_convex(a, ai, na) && _contiguous_convex(b, bj, nb))
    {
        _padic_radix_minplus_conv_lower(res, a, alen, b, blen, klo, khi);
        goto cleanup;
    }

    lambda = _best_trend(&amin, &amax, &bmin, &bmax, a, alen, b, blen);
    R = FLINT_MAX(amax - amin, bmax - bmin);

    for (i = 0; i < alen; i++)
        ad[i] = (a[i] < PADIC_RADIX_PACK_NONE) ? a[i] - lambda * i : INF;
    for (i = 0; i < blen; i++)
        bd[i] = (b[i] < PADIC_RADIX_PACK_NONE) ? b[i] - lambda * i : INF;

    if (_ks_cost(alen, blen, nout, R) <= 0.5 * _direct_cost(na, nb))
    {
        _minplus_conv_ks(res, ad, alen, amin, bd, blen, bmin, R, klo, khi);
    }
    else if (exact || _direct_cost(na, nb) <= _ks_cost(alen, blen, nout, MINPLUS_KS_MAX_RANGE))
    {
        /* (without exact, still when the bound would cost as much) */
        if ((double) na * (double) nb <= MINPLUS_BAND_PAIRS)
        {
            if (na == alen && nb == blen)
                _padic_radix_minplus_conv_exact(res, ad, alen, bd, blen, klo, khi);
            else
            {
                for (k = 0; k < nout; k++)
                    res[k] = INF;
                _minplus_sparse_update(res, ad, ai, na, bd, bj, nb, klo, khi);
            }
        }
        else
            _minplus_conv_bands(res, ad, alen, amin, amax, bd, blen, bmin, bmax, klo, khi);
    }
    else
    {
        slong delta, * q;

        result = 0;

        _padic_radix_minplus_conv_lower(res, ad, alen, bd, blen, klo, khi);

        delta = (R + MINPLUS_KS_MAX_RANGE - 1) / MINPLUS_KS_MAX_RANGE;
        q = TMP_ALLOC(sizeof(slong) * (alen + blen + nout));
        for (i = 0; i < alen; i++)
            q[i] = (ad[i] < PADIC_RADIX_PACK_NONE) ? (ad[i] - amin) / delta : INF;
        for (i = 0; i < blen; i++)
            q[alen + i] = (bd[i] < PADIC_RADIX_PACK_NONE) ? (bd[i] - bmin) / delta : INF;
        _minplus_conv_ks(q + alen + blen, q, alen, 0, q + alen, blen, 0,
            FLINT_MAX((amax - amin) / delta, (bmax - bmin) / delta), klo, khi);

        for (k = 0; k < nout; k++)
        {
            slong qk = q[alen + blen + k];
            if (qk >= PADIC_RADIX_PACK_NONE)
                res[k] = INF;       /* no finite term */
            else if (res[k] < PADIC_RADIX_PACK_NONE)
                res[k] = FLINT_MAX(res[k], amin + bmin + delta * qk);
        }

        /* outputs near the ends have few terms, and there the hull bound is
           weakest: compute those exactly */
        {
            slong T = MINPLUS_EDGE_TERMS, kend = alen + blen - 1;
            slong k0 = FLINT_MIN(khi, klo + FLINT_MAX(0, T - klo));
            slong k1 = FLINT_MAX(k0, FLINT_MAX(klo, kend - T));

            if (k0 > klo)
                _padic_radix_minplus_conv_exact(res, ad, alen, bd, blen, klo, k0);
            if (khi > k1)
                _padic_radix_minplus_conv_exact(res + (k1 - klo), ad, alen, bd, blen, k1, khi);
        }
    }

    for (k = 0; k < nout; k++)
        if (res[k] < PADIC_RADIX_PACK_NONE)
            res[k] += lambda * (k + klo);

cleanup:
    TMP_END;
    return result;
}
