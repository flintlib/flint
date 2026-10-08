/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "flint.h"
#include "ulong_extras.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "mp_real.h"
#include "impl.h"

/* Generic binary splitting for hypergeometric series in the format of
   y-cruncher's SeriesHypergeometric (its CommonP2B3 framework):

       S = (coefQ + coefP sum_{k>=1} P(k)/Q(k) prod_{j=1}^{k-1} R(j)/Q(j))
           / coefD,

   raised to the power 1 or -1, with integer polynomials P, Q, R
   (Q(k) != 0 for k >= 1, geometric convergence) and integer
   coefficients of any size.

   RECURSION.  Over a range [a, b),

       T(a,b) = sum_{k=a}^{b-1} P(k) prod_{j=a}^{k-1} R(j) prod_{j=k+1}^{b-1} Q(j),
       Q(a,b) = prod Q(k),   R(a,b) = prod R(k),

   so that the partial sum over [a, b) is T/Q, and

       T = T1 Q2 + R1 T2,   Q = Q1 Q2,   R = R1 R2,

   R being needed only by left children (the right spine skips it).
   The leaves hold L terms (about HYP_LEAF_LIMBS limbs of T, between
   HYP_LEAF_TERMS_MIN and HYP_LEAF_TERMS_MAX terms); the tree splits at
   leaf boundaries, so the term count is not rounded (except in content
   mode, below).  The
   final quotient is a single mp_real division, with coefP, coefQ, coefD
   applied as word multiplications when they fit a word.

   LEAVES.  Blocks of L terms are evaluated exactly in mpn arithmetic
   in one scratch buffer, backward in k:

       T <- P(k) Qs + R(k) T,    Qs <- Q(k) Qs,    Rs <- R(k) Rs,

   or, when R divides P (P = A R, as for pi, log 2 and most formulas),

       T <- R(k) (T + A(k) Qs),

   where A(k) is typically a word even when P(k) is not.  The
   polynomial values are computed by Horner's rule in W-word two's
   complement, with W from a bound on the values up to the final term:
   single-word values for a whole leaf up front (vectorizable), then
   the terms as mpn_mul_1 / mpn_addmul_1 passes (the sign folded into
   fused mpn_addmul_1 / mpn_submul_1, additions and copies for values
   1), in a leaf specialized to single-word values (leaf_w1); two-word
   values use umul_ppmm Horner steps and two passes, longer ones
   generic mpn code.

   CONTENT POWERS.  When Q (or R) has a large constant content c (the
   gcd of its coefficients: q^2 for atan(p/q), p^2 for p > 1), the
   products are carried without it, Q = c^(b-a) Q', and the pure powers
   c^(L 2^i) that the merges need, T = T1 Q2' c^(L 2^i) + R1 T2, come
   from a table built by squarings over a perfectly balanced tree of
   2^J leaves (the term count rounded up to a multiple of 2^J).  This
   replaces a full-size product per node (Q1 Q2) by a short one
   (Q1' Q2') but adds a full-size product by the power, and the product
   T1 Q2' does not become cheap in proportion: with FFT multiplication
   a long-by-short product costs about as much as transforming the long
   operand (measured: half of T1 Q2 for Q2' a fourteenth of Q2).  The
   net saving is a few percent at best: the table is used when the
   content has at least eight times the bits per term of the rest (at
   ratios 4-5 it was measured 3-10% slower, at 11-14 1-4% faster; never
   for the Chudnovsky constant, whose content is about half of Q), and
   the same rule decides separately for R.

   CONTENT REMOVAL.  A node only represents the ratios T/Q and R/Q, and
   g = gcd(R1, Q2) divides all of T = T1 Q2 + R1 T2, Q = Q1 Q2 and
   R = R1 R2, so R1 and Q2 can be divided by g before the products, as
   in Cheng, Hanrot, Thome, Zima and Zimmermann, "Time- and
   space-efficient evaluation of some hypergeometric constants" (ISSAC
   2007).  The cancellation is large: for zeta(3) the final Q shrinks
   from 230 to 50 bits per term, for Catalan's constant from 64 to 19,
   for log 2 from 34 to 14, for pi from 91 to 58.  The gcd is not
   computed by mpn_gcd but read off factor lists, partial
   factorizations of Q and R over the primes below 2^32 kept with every
   node: Q and R are factored over Z, the values of their linear factors
   are factored by a segmented sieve (see hyp_sieve_info) in O(sqrt(N))
   memory, as Bellard did for his pi record, and the lists of a node are
   merged from its children's.  Three things keep this cheap.  (1) A
   prime shared by R(j) (left) and Q(k) (right) divides
   u_Q u_R (k - j) + (u_R v_Q - u_Q v_R) for some pair of linear factors
   u k + v, so only the primes below a bound linear in the node's size
   are scanned, and the entries of a larger prime (most of the sieve's
   cofactors) wait in per-level buffers until the first level at which
   they can cancel instead of riding up through every merge.  (2) The
   removal is done where it pays (_hyp_remove_pays): an exact division
   costs up to about two products of its size, against the products it
   shrinks at this node and every node above; near the root it does not
   pay, and the lists are dropped from there up.  (3) Only exact values
   take part (the top levels are truncated to the working precision).
   At smaller precisions the sieve and the lists cost about as much as
   the removal saves, and it is not done (_hyp_gcd_min_limbs).

   TERMS AND TAIL.  The tail after N terms is bounded rigorously (see
   tail_struct): the product of the ratios |R(j)/Q(j)| is accumulated
   exactly (up to double rounding, bounded) for small j and bounded in
   closed form beyond, and the terms beyond N by a geometric series
   from leading-term bounds.  The term count is the least N whose bound
   is below the target, and the bound is added to the ball. */

#define HYP_LEAF_LIMBS 24
#define HYP_LEAF_TERMS_MAX 32
/* at least this many terms per leaf even when the terms are large (as for
   Zuniga series with multi-limb coefficients): short leaves of long terms
   cost more in merges than they save in the leaf */
#define HYP_LEAF_TERMS_MIN 24
/* the same with the content removal */
#define HYP_GCD_LEAF_LIMBS 40
#define HYP_GCD_LEAF_TERMS_MIN 12

/* ---- polynomials ---- */

/* the values of an integer polynomial at 1 <= k <= kmax, by Horner's
   rule in W-word two's complement, exact since every Horner
   intermediate is below 2^(FLINT_BITS W - 1) in absolute value */
typedef struct
{
    slong len;          /* degree + 1 (0: the zero polynomial) */
    slong W;
    nn_ptr c;           /* len * W limbs, two's complement */
}
hpoly_struct;

typedef hpoly_struct hpoly_t[1];

/* an upper bound for log2 sum_i |f_i| k^i (from the bit sizes, for
   coefficients of any size) */
static double
_log2_bound_fmpz(const fmpz * f, slong len, ulong k)
{
    double m = 0.0, lk = log2((double) FLINT_MAX(k, 1));
    slong i;

    for (i = 0; i < len; i++)
        if (!fmpz_is_zero(f + i))
            m = FLINT_MAX(m, (double) fmpz_bits(f + i) + i * lk);
    return (m + FLINT_BIT_COUNT(len)) * (1.0 + 1e-12) + 1e-6;
}

static void
hpoly_init(hpoly_t p, const fmpz * f, slong len, ulong kmax)
{
    while (len > 0 && fmpz_is_zero(f + len - 1))
        len--;

    p->len = len;

    /* log2 of sum |c_i| kmax^i, an upper bound for every Horner
       intermediate */
    p->W = (slong) ((_log2_bound_fmpz(f, len, kmax) + 2.0) / FLINT_BITS) + 1;
    p->c = NULL;
}

/* store the coefficients in c (room for len W limbs) */
static void
hpoly_fill(hpoly_t p, const fmpz * f, nn_ptr c)
{
    slong i, j, len = p->len;

    p->c = c;
    for (i = 0; i < len; i++)
    {
        nn_ptr d = p->c + i * p->W;

        for (j = 0; j < p->W; j++)
            d[j] = 0;
        if (!COEFF_IS_MPZ(f[i]))
        {
            d[0] = FLINT_UABS(f[i]);
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_abs(t, f + i);
            FLINT_ASSERT(fmpz_size(t) <= p->W);
            fmpz_get_ui_array(d, p->W, t);
            fmpz_clear(t);
        }
        if (fmpz_sgn(f + i) < 0)
            mpn_neg(d, d, p->W);
    }
}

/* (v, *vn) = |p(k)|, returning the sign (1 for negative); v has room for
   W limbs; *vn = 0 for p(k) = 0 */
static int
hpoly_eval(nn_ptr v, slong * vn, const hpoly_t p, ulong k)
{
    slong i, W = p->W, len = p->len;
    int neg;

    if (len == 0)
    {
        *vn = 0;
        return 0;
    }

    if (W == 1)
    {
        slong x = (slong) p->c[len - 1];
        for (i = len - 2; i >= 0; i--)
            x = x * (slong) k + (slong) p->c[i];
        neg = (x < 0);
        v[0] = FLINT_UABS(x);
        *vn = (v[0] != 0);
        return neg;
    }
    else if (W == 2)
    {
        ulong hi = p->c[2 * (len - 1) + 1], lo = p->c[2 * (len - 1)];
        for (i = len - 2; i >= 0; i--)
        {
            ulong ph, pl;
            umul_ppmm(ph, pl, lo, k);
            ph += hi * k;
            add_ssaaaa(hi, lo, ph, pl, p->c[2 * i + 1], p->c[2 * i]);
        }
        neg = ((slong) hi < 0);
        if (neg)
            sub_ddmmss(hi, lo, UWORD(0), UWORD(0), hi, lo);
        v[0] = lo;
        v[1] = hi;
        *vn = (hi != 0) ? 2 : (lo != 0);
        return neg;
    }
    else
    {
        flint_mpn_copyi(v, p->c + (len - 1) * W, W);
        for (i = len - 2; i >= 0; i--)
        {
            mpn_mul_1(v, v, W, k);
            mpn_add_n(v, v, p->c + i * W, W);
        }
        neg = ((slong) v[W - 1] < 0);
        if (neg)
            mpn_neg(v, v, W);
        i = W;
        while (i > 0 && v[i - 1] == 0)
            i--;
        *vn = i;
        return neg;
    }
}

/* the values of p (W = 1) at k = a, ..., a + m - 1, by Horner's rule
   in wrapping word arithmetic (exact: every intermediate fits) */
static void
hpoly_eval_block(ulong * v, const hpoly_t p, slong a, slong m)
{
    slong i, j, len = p->len;

    if (len == 0)
    {
        for (j = 0; j < m; j++)
            v[j] = 0;
        return;
    }

    for (j = 0; j < m; j++)
        v[j] = p->c[len - 1];
    for (i = len - 2; i >= 0; i--)
    {
        ulong c = p->c[i];
        for (j = 0; j < m; j++)
            v[j] = v[j] * (ulong) (a + j) + c;
    }
}

/* (v, *vn) = |p(k)| with the sign returned, from the block values when
   W = 1 */
static inline int
hpoly_value(nn_ptr v, slong * vn, const hpoly_t p, const ulong * blk,
    slong i, slong k)
{
    if (p->W == 1)
    {
        slong x = (slong) blk[i];
        v[0] = FLINT_UABS(x);
        *vn = (x != 0);
        return x < 0;
    }
    return hpoly_eval(v, vn, p, k);
}

/* ---- leaf arithmetic on (limbs, length, sign) ---- */

static inline slong
_normalize(nn_srcptr x, slong xn)
{
    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    return xn;
}

/* z = x v, v of vn >= 1 limbs; z must not alias x unless vn = 1;
   returns the length */
static inline slong
_mul_small(nn_ptr z, nn_srcptr x, slong xn, nn_srcptr v, slong vn)
{
    ulong cy;

    if (xn == 0)
        return 0;
    if (vn == 1)
    {
        cy = mpn_mul_1(z, x, xn, v[0]);
        z[xn] = cy;
        return xn + (cy != 0);
    }
    if (vn == 2)
    {
        z[xn] = mpn_mul_1(z, x, xn, v[0]);
        z[xn + 1] = mpn_addmul_1(z + 1, x, xn, v[1]);
        return _normalize(z, xn + 2);
    }
    if (xn >= vn)
        flint_mpn_mul(z, x, xn, v, vn);
    else
        flint_mpn_mul(z, v, vn, x, xn);
    return _normalize(z, xn + vn);
}

/* (x, *xn, *xneg) += (-1)^yneg (y, yn) w for a single limb w; x has
   room for max(xn, yn) + 1 limbs */
static inline void
_addmul_1_signed(nn_ptr x, slong * xn, int * xneg, nn_srcptr y, slong yn,
    ulong w, int yneg)
{
    slong n = *xn;
    ulong cy;

    if (yn == 0 || w == 0)
        return;

    if (n == 0)
    {
        cy = mpn_mul_1(x, y, yn, w);
        x[yn] = cy;
        *xn = yn + (cy != 0);
        *xneg = yneg;
        return;
    }

    /* opposite signs with |y| > |x| and w = 1: y - x directly */
    if (w == 1 && *xneg != yneg
        && (yn > n || (yn == n && mpn_cmp(y, x, n) > 0)))
    {
        mpn_sub(x, y, yn, x, n);
        *xn = _normalize(x, yn);
        *xneg = yneg;
        return;
    }

    if (n < yn)
    {
        flint_mpn_zero(x + n, yn - n);
        n = yn;
    }

    if (*xneg == yneg)
    {
        cy = (w == 1) ? mpn_add_n(x, x, y, yn) : mpn_addmul_1(x, y, yn, w);
        if (n > yn)
            cy = mpn_add_1(x + yn, x + yn, n - yn, cy);
        x[n] = cy;
        *xn = n + (cy != 0);
    }
    else
    {
        cy = (w == 1) ? mpn_sub_n(x, x, y, yn) : mpn_submul_1(x, y, yn, w);
        if (n > yn)
            cy = mpn_sub_1(x + yn, x + yn, n - yn, cy);
        if (cy != 0)
        {
            /* the stored value minus cy B^n is negative: its absolute
               value is (cy - [x != 0]) B^n + (B^n - x) */
            cy -= mpn_neg(x, x, n);
            x[n] = cy;
            *xneg = !*xneg;
            n += (cy != 0);
        }
        *xn = _normalize(x, n);
    }
}

/* (x, *xn, *xneg) += (-1)^yneg (y, yn); x has room for
   max(xn, yn) + 1 limbs */
static void
_add_signed(nn_ptr x, slong * xn, int * xneg, nn_srcptr y, slong yn,
    int yneg)
{
    slong L;

    if (yn == 0)
        return;

    L = FLINT_MAX(*xn, yn);
    if (*xn < L)
        flint_mpn_zero(x + *xn, L - *xn);

    if (*xn == 0)
        *xneg = yneg;

    if (*xneg == yneg)
    {
        x[L] = mpn_add(x, x, L, y, yn);
        *xn = _normalize(x, L + 1);
    }
    else
    {
        if (mpn_sub(x, x, L, y, yn))
        {
            mpn_neg(x, x, L);
            *xneg = !*xneg;
        }
        *xn = _normalize(x, L);
    }
}

/* ---- content removal: factor lists and the sieve ---- */

/* A factor list is a partial factorization of an exact integer X: a list
   of (p, e) sorted by p, the product of the p^e dividing X, with p < 2^32
   a prime or an unsplit cofactor (see the sieve; either way the entries
   come from disjoint parts of X), and exponents saturating at
   2^32 - 1, so that sums and differences keep them lower bounds.  An
   integer whose list is not tracked (ok = 0) takes part in no
   removal. */
typedef struct
{
    uint32_t p, e;
}
hyp_fe;

typedef struct
{
    hyp_fe * d;
    slong len, alloc;
    int ok;
}
hyp_flist;

static void
flist_init(hyp_flist * l)
{
    l->d = NULL;
    l->len = l->alloc = 0;
    l->ok = 0;
}

static void
flist_clear(hyp_flist * l)
{
    flint_free(l->d);
    flist_init(l);
}

static void
flist_fit(hyp_flist * l, slong n)
{
    if (l->alloc < n)
    {
        n = FLINT_MAX(n, 2 * l->alloc);
        l->d = flint_realloc(l->d, n * sizeof(hyp_fe));
        l->alloc = n;
    }
}

static inline uint32_t
_esat(uint64_t e)
{
    return (e > UINT32_MAX) ? UINT32_MAX : (uint32_t) e;
}

/* a <- a + b (exponents of the product), entries of exponent 0 (left
   by flist_gcd_remove) dropped; b unchanged; the result is formed in
   the scratch list t and swapped into a */
static void
flist_add(hyp_flist * a, const hyp_flist * b, hyp_flist * t)
{
    slong i, j, n, na, nb;
    const hyp_fe * A, * Bp;
    hyp_fe * r;

    if (!a->ok || !b->ok)
    {
        a->ok = 0;
        a->len = 0;
        return;
    }

    na = a->len;
    nb = b->len;
    flist_fit(t, na + nb + 1);
    A = a->d;
    Bp = b->d;
    r = t->d;

    i = j = n = 0;
    if (nb == 0)
    {
        for (i = 0; i < na; i++)
        {
            r[n] = A[i];
            n += (A[i].e != 0);
        }
    }
    else
    {
        /* branch-free (the merge of random primes defeats prediction):
           both heads loaded, the exponents masked */
        while (i < na && j < nb)
        {
            hyp_fe x = A[i], y = Bp[j];
            uint32_t p = FLINT_MIN(x.p, y.p);
            uint64_t ta = (x.p == p), tb = (y.p == p);
            uint64_t e = ((uint64_t) x.e & (0 - ta))
                + ((uint64_t) y.e & (0 - tb));
            r[n].p = p;
            r[n].e = (uint32_t) FLINT_MIN(e, (uint64_t) UINT32_MAX);
            n += (e != 0);
            i += ta;
            j += tb;
        }
        for ( ; i < na; i++)
        {
            r[n] = A[i];
            n += (A[i].e != 0);
        }
        for ( ; j < nb; j++)
        {
            r[n] = Bp[j];
            n += (Bp[j].e != 0);
        }
    }

    FLINT_SWAP(hyp_fe *, a->d, t->d);
    FLINT_SWAP(slong, a->alloc, t->alloc);
    a->len = n;
}

/* a <- a + b + c with c short (merged into the pass as a third
   source), entries of exponent 0 dropped; the result is formed in t */
static void
flist_add3(hyp_flist * a, const hyp_flist * b, const hyp_flist * c,
    hyp_flist * t)
{
    slong i = 0, j = 0, k = 0, n = 0, na, nb, nc;
    const hyp_fe * A, * Bp, * C;
    hyp_fe * r;

    if (c->len == 0)
    {
        flist_add(a, b, t);
        return;
    }
    if (!a->ok || !b->ok)
    {
        a->ok = 0;
        a->len = 0;
        return;
    }

    na = a->len;
    nb = b->len;
    nc = c->len;
    flist_fit(t, na + nb + nc);
    A = a->d;
    Bp = b->d;
    C = c->d;
    r = t->d;

    /* the two long lists by the two-way loop; the short one is checked
       against the smaller head, a branch rarely taken */
    {
        uint32_t pc = C[0].p;
        while (i < na && j < nb)
        {
            hyp_fe x = A[i], y = Bp[j];
            uint32_t p = FLINT_MIN(x.p, y.p);
            uint64_t ta = (x.p == p), tb = (y.p == p), e;
            if (pc <= p)
            {
                if (pc < p)
                {
                    r[n++] = C[k++];
                    pc = (k < nc) ? C[k].p : UINT32_MAX;
                    continue;
                }
                e = C[k++].e;
                pc = (k < nc) ? C[k].p : UINT32_MAX;
            }
            else
                e = 0;
            e += ((uint64_t) x.e & (0 - ta)) + ((uint64_t) y.e & (0 - tb));
            r[n].p = p;
            r[n].e = (uint32_t) FLINT_MIN(e, (uint64_t) UINT32_MAX);
            n += (e != 0);
            i += ta;
            j += tb;
        }
        /* the rest, three-way */
        while (i < na || j < nb || k < nc)
        {
            uint32_t pa = (i < na) ? A[i].p : UINT32_MAX;
            uint32_t pb = (j < nb) ? Bp[j].p : UINT32_MAX;
            uint32_t pcc = (k < nc) ? C[k].p : UINT32_MAX;
            uint32_t p = FLINT_MIN(FLINT_MIN(pa, pb), pcc);
            uint64_t e = 0;
            if (pa == p)
                e += A[i++].e;
            if (pb == p)
                e += Bp[j++].e;
            if (pcc == p)
                e += C[k++].e;
            r[n].p = p;
            r[n].e = _esat(e);
            n += (e != 0);
        }
    }

    FLINT_SWAP(hyp_fe *, a->d, t->d);
    FLINT_SWAP(slong, a->alloc, t->alloc);
    a->len = n;
}

/* g = gcd(a + a2, b + b2) entrywise over the primes p <= bound, removed
   from the sources (a before a2, b before b2), leaving entries of
   exponent 0 (dropped by the next flist_add) */
static void
flist_gcd_remove(hyp_flist * g, hyp_flist * a, hyp_flist * a2,
    hyp_flist * b, hyp_flist * b2, ulong bound)
{
    slong i = 0, i2 = 0, j = 0, j2 = 0;
    slong na = a->len, na2 = a2->len, nb = b->len, nb2 = b2->len;
    hyp_fe * A = a->d, * A2 = a2->d, * Bp = b->d, * B2 = b2->d;

    g->len = 0;
    if (!a->ok || !b->ok)
        return;

    for (;;)
    {
        uint32_t pa = (i < na) ? A[i].p : UINT32_MAX;
        uint32_t qa = (i2 < na2) ? A2[i2].p : UINT32_MAX;
        uint32_t pb = (j < nb) ? Bp[j].p : UINT32_MAX;
        uint32_t qb = (j2 < nb2) ? B2[j2].p : UINT32_MAX;
        uint32_t ra = FLINT_MIN(pa, qa), rb = FLINT_MIN(pb, qb);

        /* UINT32_MAX (not a prime) marks an exhausted side */
        if (ra > bound || rb > bound || ra == UINT32_MAX
            || rb == UINT32_MAX)
            break;
        if (ra < rb)
        {
            i += (pa == ra);
            i2 += (qa == ra);
        }
        else if (rb < ra)
        {
            j += (pb == rb);
            j2 += (qb == rb);
        }
        else
        {
            uint64_t ea = ((pa == ra) ? A[i].e : 0)
                + (uint64_t) ((qa == ra) ? A2[i2].e : 0);
            uint64_t eb = ((pb == rb) ? Bp[j].e : 0)
                + (uint64_t) ((qb == rb) ? B2[j2].e : 0);
            uint64_t m = FLINT_MIN(ea, eb), x;

            if (m != 0)
            {
                flist_fit(g, g->len + 1);
                g->d[g->len].p = ra;
                g->d[g->len++].e = _esat(m);
                /* subtract m: a first, then a2; b first, then b2 */
                x = m;
                if (pa == ra)
                {
                    uint32_t y = (uint32_t) FLINT_MIN(x, (uint64_t) A[i].e);
                    A[i].e -= y;
                    x -= y;
                }
                if (x != 0 && qa == ra)
                    A2[i2].e -= (uint32_t) FLINT_MIN(x, (uint64_t) A2[i2].e);
                x = m;
                if (pb == rb)
                {
                    uint32_t y = (uint32_t) FLINT_MIN(x, (uint64_t) Bp[j].e);
                    Bp[j].e -= y;
                    x -= y;
                }
                if (x != 0 && qb == rb)
                    B2[j2].e -= (uint32_t) FLINT_MIN(x, (uint64_t) B2[j2].e);
            }
            i += (pa == ra);
            i2 += (qa == ra);
            j += (pb == rb);
            j2 += (qb == rb);
        }
    }
}

/* the product of n words, n >= 1, by a balanced tree; (z, return value)
   with z room for n limbs; t scratch of n limbs */
static slong
_words_prod(nn_ptr z, const ulong * w, slong n, nn_ptr t)
{
    slong zn, m, an, bn;

    if (n <= 16)
    {
        z[0] = w[0];
        zn = 1;
        for (m = 1; m < n; m++)
        {
            ulong cy = mpn_mul_1(z, z, zn, w[m]);
            z[zn] = cy;
            zn += (cy != 0);
        }
        return zn;
    }

    m = n / 2;
    an = _words_prod(t, w, m, z);
    bn = _words_prod(t + m, w + m, n - m, z);
    if (an >= bn)
        flint_mpn_mul(z, t, an, t + m, bn);
    else
        flint_mpn_mul(z, t + m, bn, t, an);
    zn = an + bn;
    return zn - (z[zn - 1] == 0);
}

/* an upper estimate of log2 prod p^e over the list */
static double
flist_bits(const hyp_flist * g)
{
    double b = 0.0;
    slong i;
    for (i = 0; i < g->len; i++)
        b += (double) g->d[i].e * FLINT_BIT_COUNT(g->d[i].p);
    return b;
}

#define HYP_SMALL_G 64

/* prod p^e over the list, as (G, return value): in sbuf (HYP_SMALL_G
   limbs) when it fits there, else freshly allocated */
static slong
flist_expand(nn_ptr * G, const hyp_flist * g, nn_ptr sbuf)
{
    slong i, n = 0, alloc = HYP_SMALL_G, gn;
    ulong wbuf[HYP_SMALL_G], * w = wbuf, acc = 1, hi, lo, k;
    nn_ptr t;

    for (i = 0; i < g->len; i++)
    {
        ulong p = g->d[i].p;
        for (k = 0; k < g->d[i].e; k++)
        {
            umul_ppmm(hi, lo, acc, p);
            if (hi == 0)
                acc = lo;
            else
            {
                if (n == alloc)
                {
                    alloc *= 2;
                    if (w == wbuf)
                    {
                        w = flint_malloc(alloc * sizeof(ulong));
                        memcpy(w, wbuf, n * sizeof(ulong));
                    }
                    else
                        w = flint_realloc(w, alloc * sizeof(ulong));
                }
                w[n++] = acc;
                acc = p;
            }
        }
    }
    if (n == alloc)
    {
        if (w == wbuf)
        {
            w = flint_malloc((alloc + 1) * sizeof(ulong));
            memcpy(w, wbuf, n * sizeof(ulong));
        }
        else
            w = flint_realloc(w, (alloc + 1) * sizeof(ulong));
    }
    w[n++] = acc;

    if (n <= HYP_SMALL_G / 2)
    {
        *G = sbuf;
        gn = _words_prod(sbuf, w, n, sbuf + n);
    }
    else
    {
        *G = flint_malloc(n * sizeof(ulong));
        t = flint_malloc((n + 1) * sizeof(ulong));
        gn = _words_prod(*G, w, n, t);
        flint_free(t);
    }
    if (w != wbuf)
        flint_free(w);
    return gn;
}

/* x / G for an exact integer x divisible by G, with the precomputed
   inverse pre when not NULL */
static void
_mp_real_divexact_mpn_pre(mp_real_t x, nn_srcptr G, slong gn,
    const flint_mpn_divexact_preinv_struct * pre)
{
    slong z = x->exp - x->size, an = x->exp, qn;
    nn_ptr a, q;
    int neg = x->negative;
    ulong sbuf[2 * HYP_SMALL_G + 2];

    if (gn == 1 && G[0] == 1)
        return;

    FLINT_ASSERT(x->err == 0 && z >= 0 && an >= gn);

    qn = an - gn + 1;
    if (qn + (z ? an : 0) <= 2 * HYP_SMALL_G + 2)
        q = sbuf;
    else
        q = flint_malloc((qn + (z ? an : 0)) * sizeof(ulong));
    if (z == 0)
        a = x->d;
    else
    {
        a = q + qn;
        flint_mpn_zero(a, z);
        flint_mpn_copyi(a + z, x->d, x->size);
    }
    if (pre != NULL)
        flint_mpn_divexact_preinv(q, a, an, pre);
    else
        flint_mpn_divexact(q, a, an, G, gn);
    _mp_real_set_mpn_2exp(x, q, qn, 0);
    if (neg)
        mp_real_neg(x, x);
    if (q != sbuf)
        flint_free(q);
}

static inline int
_mp_real_is_exact_int(const mp_real_t x)
{
    return x->size != 0 && x->err == 0 && x->exp >= x->size;
}

/* The sieve.  Q and R are factored over Z; their linear factors
   f_i(k) = u_i k + v_i (u_i > 0, gcd(u_i, v_i) = 1, of multiplicity
   mQ_i in Q and mR_i in R) are sieved over windows of consecutive
   terms with the primes p <= sqrt(M), M the largest |f_i(k)|, from
   the roots -v_i / u_i mod p, and the cofactors left (primes beyond
   sqrt(M)) are recorded too; the contents are factored by trial
   division (a cofactor left unfactored is just not tracked, as are
   factors of degree > 1).  The memory is O(sqrt(M)): the primes and
   roots (read-only, shared by the threads) and one window per thread,
   refilled as the depth-first traversal reaches its leaves in
   increasing order, of a few thousand terms (enough that the roots,
   visited once per window, cost a few operations per term) within
   about 8 MB; there is no table of all the values as in the sieve of
   gmp-chudnovsky. */
typedef struct
{
    slong nf;
    ulong * u;
    slong * v;
    uint32_t * mQ, * mR;
    slong np;
    uint32_t * primes;
    ulong * pinv, * plim;       /* p^-1 mod 2^FLINT_BITS, floor(UWORD_MAX/p) */
    uint32_t * roots;           /* np * nf; UINT32_MAX: p | u_i */
    hyp_flist cq, cr;           /* content primes */
    slong L, N, nleaves, segleaves;
    /* a prime shared by R(j) and Q(k), j < k, from linear factors
       divides u_Q u_R (k - j) + (u_R v_Q - u_Q v_R), so it is at most
       cmul (k - j) + cadd unless that vanishes (a telescoping pair:
       nobound), or it is a content prime (at most cmaxp) */
    double cmul, cadd, cmaxp;
    int nobound;
    double lb[64];              /* the bound at a node of 2^lev leaves,
                                   lev < 63 (see HYP_LB_INV) */
    slong cQ, cR, dQ, dR;       /* per-leaf capacities of a window */
    /* the same bound for the entries of one factor, as an entry of Q
       (lbq) or of R (lbr): nf rows of 64, like lb */
    double * lbq, * lbr;
    /* an entry of Q (R) with a prime beyond every R (Q) value cannot
       cancel */
    ulong qprune, rprune;
}
hyp_sieve_info;

/* one thread's window: leaves j0 <= j < j1.  For the leaf of index l =
   j - j0, the sieved primes form a list sorted by p, sQ + l cQ (nQ[l]
   entries): the sieve visits the primes in increasing order and appends
   to the leaf's list, or adds to its last entry; the cofactors (primes
   beyond the sieve bound, nearly all deferred to higher levels) are
   kept unsorted in fQ + l dQ (mQ[l] entries).  The same for R. */
typedef struct
{
    slong j0, j1;
    ulong * vals;               /* nf * (terms of the window) */
    hyp_fe * sQ, * sR, * fQ, * fR;
    uint32_t * nQ, * nR, * mQ, * mR;
}
hyp_sieve_seg;

#define HYP_NO_ROOT UINT32_MAX

static void
seg_init(hyp_sieve_seg * S)
{
    memset(S, 0, sizeof(hyp_sieve_seg));
}

static void
seg_clear(hyp_sieve_seg * S)
{
    flint_free(S->vals);
    flint_free(S->sQ);
    flint_free(S->nQ);
    seg_init(S);
}

/* sort the triples (leaf, p, e) ev[lo, hi) by p, stably: by insertion
   when short, else by LSD radix (8- or 11-bit digits up to the top bit
   of the largest p) */
static void
_ev_radix_sort_p(uint32_t * ev, slong lo, slong hi)
{
    slong n = hi - lo, i, d, nd;
    uint32_t maxp = 0, * a = ev + 3 * lo, * b, * t, mask;
    slong cnt[2048];
    int shift, bits;

    if (n < 2)
        return;

    /* short: insertion sort (stable) */
    if (n <= 48)
    {
        for (i = 1; i < n; i++)
        {
            uint32_t x0 = a[3 * i], x1 = a[3 * i + 1], x2 = a[3 * i + 2];
            slong j = i;
            while (j > 0 && a[3 * (j - 1) + 1] > x1)
            {
                a[3 * j] = a[3 * (j - 1)];
                a[3 * j + 1] = a[3 * (j - 1) + 1];
                a[3 * j + 2] = a[3 * (j - 1) + 2];
                j--;
            }
            a[3 * j] = x0;
            a[3 * j + 1] = x1;
            a[3 * j + 2] = x2;
        }
        return;
    }

    for (i = 0; i < n; i++)
        maxp = FLINT_MAX(maxp, a[3 * i + 1]);

    bits = (n < 4096) ? 8 : 11;
    nd = WORD(1) << bits;
    mask = (uint32_t) nd - 1;
    b = flint_malloc(3 * n * sizeof(uint32_t));
    t = a;
    for (shift = 0; shift < 32 && (maxp >> shift) != 0; shift += bits)
    {
        slong pos = 0;
        memset(cnt, 0, nd * sizeof(slong));
        for (i = 0; i < n; i++)
            cnt[(t[3 * i + 1] >> shift) & mask]++;
        for (d = 0; d < nd; d++)
        {
            slong c = cnt[d];
            cnt[d] = pos;
            pos += c;
        }
        for (i = 0; i < n; i++)
        {
            slong j = cnt[(t[3 * i + 1] >> shift) & mask]++;
            b[3 * j] = t[3 * i];
            b[3 * j + 1] = t[3 * i + 1];
            b[3 * j + 2] = t[3 * i + 2];
        }
        FLINT_SWAP(uint32_t *, t, b);
    }
    if (t != a)
    {
        memcpy(a, t, 3 * n * sizeof(uint32_t));
        b = t;
    }
    flint_free(b);
}

/* append (p, e) to a leaf's list of capacity cap, or add to its last
   entry; an entry beyond the capacity is dropped (a list is a lower
   bound, so this only loses some cancellation) */
static inline void
_leaf_append(hyp_fe * lst, uint32_t * n, uint32_t p, uint64_t e, slong cap)
{
    uint32_t k = *n;
    if (k != 0 && lst[k - 1].p == p)
        lst[k - 1].e = _esat((uint64_t) lst[k - 1].e + e);
    else if (k < cap)
    {
        lst[k].p = p;
        lst[k].e = _esat(e);
        *n = k + 1;
    }
}

static inline slong _hyp_level_row(const hyp_sieve_info * I,
    const double * lb, uint32_t p);

/* sieve the window of the leaves j0, ... (up to segleaves of them) */
static void
seg_fill(hyp_sieve_seg * S, const hyp_sieve_info * I, slong j0)
{
    slong j1 = FLINT_MIN(j0 + I->segleaves, I->nleaves);
    slong L = I->L, nf = I->nf, nl = j1 - j0, i, t, pi, l, c;
    slong k0 = 1 + j0 * L, k1 = FLINT_MIN(1 + j1 * L, I->N + 1), K = k1 - k0;
    slong cQ = I->cQ, cR = I->cR, dQ = I->dQ, dR = I->dR;
    slong cqi = 0, cri = 0;

    if (S->vals == NULL)
    {
        slong sl = I->segleaves;
        S->vals = flint_malloc((nf * sl * L + 1) * sizeof(ulong));
        S->sQ = flint_malloc((sl * (cQ + cR + dQ + dR) + 1) * sizeof(hyp_fe));
        S->sR = S->sQ + sl * cQ;
        S->fQ = S->sR + sl * cR;
        S->fR = S->fQ + sl * dQ;
        S->nQ = flint_malloc((4 * sl + 1) * sizeof(uint32_t));
        S->nR = S->nQ + sl;
        S->mQ = S->nR + sl;
        S->mR = S->mQ + sl;
    }

    S->j0 = j0;
    S->j1 = j1;
    memset(S->nQ, 0, 4 * I->segleaves * sizeof(uint32_t));

    for (i = 0; i < nf; i++)
    {
        ulong * vi = S->vals + i * K;
        for (t = 0; t < K; t++)
        {
            slong x = (slong) I->u[i] * (k0 + t) + I->v[i];
            vi[t] = FLINT_UABS(x);
        }
    }

    for (pi = 0; pi < I->np; pi++)
    {
        ulong p = I->primes[pi], pinv = I->pinv[pi], plim = I->plim[pi];
        ulong k0p = (ulong) k0 % p, stl = p / L, sto = p % L;
        const uint32_t * rt = I->roots + pi * nf;
        int okq = (p <= I->qprune), okr = (p <= I->rprune);

        for (i = 0; i < nf; i++)
        {
            ulong * vi = S->vals + i * K;
            ulong r = rt[i];
            uint32_t mq = okq ? I->mQ[i] : 0, mr = okr ? I->mR[i] : 0;
            slong off;

            if (r == HYP_NO_ROOT)
                continue;
            t = (slong) ((r >= k0p) ? r - k0p : r + p - k0p);
            l = t / L;
            off = t % L;
            for ( ; t < K; t += p)
            {
                ulong x = vi[t];
                uint64_t e = 0;
                if (p == 2)
                {
                    e = flint_ctz(x);
                    x >>= e;
                }
                else
                {
                    do
                    {
                        x *= pinv;
                        e++;
                    }
                    while (x * pinv <= plim);
                }
                vi[t] = x;
                if (mq)
                    _leaf_append(S->sQ + l * cQ, S->nQ + l, (uint32_t) p,
                        e * mq, cQ);
                if (mr)
                    _leaf_append(S->sR + l * cR, S->nR + l, (uint32_t) p,
                        e * mr, cR);
                l += stl;
                off += sto;
                if (off >= L)
                {
                    off -= L;
                    l++;
                }
            }
        }

        /* content primes among the sieve primes */
        if (cqi < I->cq.len && I->cq.d[cqi].p == p)
        {
            for (l = 0; l < nl; l++)
                _leaf_append(S->sQ + l * cQ, S->nQ + l, (uint32_t) p,
                    (uint64_t) I->cq.d[cqi].e * FLINT_MIN(L, K - l * L), cQ);
            cqi++;
        }
        if (cri < I->cr.len && I->cr.d[cri].p == p)
        {
            for (l = 0; l < nl; l++)
                _leaf_append(S->sR + l * cR, S->nR + l, (uint32_t) p,
                    (uint64_t) I->cr.d[cri].e * FLINT_MIN(L, K - l * L), cR);
            cri++;
        }
    }

    /* content primes beyond the sieve primes */
    for (c = cqi; c < I->cq.len; c++)
        for (l = 0; l < nl; l++)
            _leaf_append(S->sQ + l * cQ, S->nQ + l, I->cq.d[c].p,
                (uint64_t) I->cq.d[c].e * FLINT_MIN(L, K - l * L), cQ);
    for (c = cri; c < I->cr.len; c++)
        for (l = 0; l < nl; l++)
            _leaf_append(S->sR + l * cR, S->nR + l, I->cr.d[c].p,
                (uint64_t) I->cr.d[c].e * FLINT_MIN(L, K - l * L), cR);

    /* the cofactors: primes beyond the sieve bound, written branch-free
       (each value is stored and kept or not; the slots have one to
       spare): whether a value has one is unpredictable */
    for (i = 0; i < nf; i++)
    {
        const ulong * vi = S->vals + i * K;
        uint32_t mq = I->mQ[i], mr = I->mR[i];
        const double * lq = I->lbq + 64 * i, * lr = I->lbr + 64 * i;
        for (l = 0, t = 0; l < nl; l++)
        {
            slong te = FLINT_MIN(t + L, K);
            if (mq)
            {
                hyp_fe * d = S->fQ + l * dQ;
                uint32_t k = S->mQ[l];
                slong u;
                for (u = t; u < te; u++)
                {
                    ulong x = vi[u];
                    d[k].p = (uint32_t) x;
                    d[k].e = mq | ((uint32_t) _hyp_level_row(I, lq,
                        (uint32_t) x) << 24);
                    k += (x > 1) & (x <= I->qprune);
                }
                S->mQ[l] = k;
            }
            if (mr)
            {
                hyp_fe * d = S->fR + l * dR;
                uint32_t k = S->mR[l];
                slong u;
                for (u = t; u < te; u++)
                {
                    ulong x = vi[u];
                    d[k].p = (uint32_t) x;
                    d[k].e = mr | ((uint32_t) _hyp_level_row(I, lr,
                        (uint32_t) x) << 24);
                    k += (x > 1) & (x <= I->rprune);
                }
                S->mR[l] = k;
            }
            t = te;
        }
    }
}

/* Deferral.  At a node of m terms, a prime of an R entry of the left
   child and a Q entry of the right child divides u_Q u_R d + w for some
   d < m (or is a content prime), so the entries of a prime p cannot
   cancel below the first level lev with cmul 2^lev L + cadd >= p.  The
   leaf lists keep the entries that can cancel at level 1 (most of the
   small primes); the others (most of the cofactors) wait in per-level
   buffers, as (leaf, p, e), and join the lists of the two children of
   the node of their level just before its removal.  Without this, the
   cofactors would be carried through every merge from the leaves up,
   which dominated the cost of the lists. */
typedef struct
{
    uint32_t * d;               /* triples (leaf, p, e) */
    slong len, alloc;
}
hyp_pend;

typedef struct
{
    hyp_sieve_seg seg;
    hyp_pend * pq, * pr;        /* per level, 0 <= lev <= nlev - 1 */
    slong nlev;
    slong lc;                   /* entries for levels >= lc are discarded */
    hyp_flist scr, ins;         /* scratch lists */
    hyp_flist pl[4];            /* a merge's pending lists */
}
hyp_thread;

static void
pend_push(hyp_pend * P, uint32_t leaf, uint32_t p, uint32_t e)
{
    if (P->len == P->alloc)
    {
        P->alloc = FLINT_MAX(16, 2 * P->alloc);
        P->d = flint_realloc(P->d, 3 * P->alloc * sizeof(uint32_t));
    }
    P->d[3 * P->len] = leaf;
    P->d[3 * P->len + 1] = p;
    P->d[3 * P->len + 2] = e;
    P->len++;
}

static void
pend_append(hyp_pend * P, const hyp_pend * Q)
{
    if (Q->len == 0)
        return;
    if (P->len + Q->len > P->alloc)
    {
        P->alloc = FLINT_MAX(P->len + Q->len, 2 * P->alloc);
        P->d = flint_realloc(P->d, 3 * P->alloc * sizeof(uint32_t));
    }
    memcpy(P->d + 3 * P->len, Q->d, 3 * Q->len * sizeof(uint32_t));
    P->len += Q->len;
}

static void
pend_free(hyp_pend * P)
{
    flint_free(P->d);
    P->d = NULL;
    P->len = P->alloc = 0;
}

static void
thread_init(hyp_thread * t, slong nlev)
{
    slong i;
    seg_init(&t->seg);
    t->nlev = nlev;
    t->lc = nlev;
    t->pq = flint_calloc(2 * nlev, sizeof(hyp_pend));
    t->pr = t->pq + nlev;
    flist_init(&t->scr);
    flist_init(&t->ins);
    for (i = 0; i < 4; i++)
        flist_init(t->pl + i);
}

static void
thread_clear(hyp_thread * t)
{
    slong i;
    seg_clear(&t->seg);
    for (i = 0; i < 2 * t->nlev; i++)
        pend_free(t->pq + i);
    flint_free(t->pq);
    flist_clear(&t->scr);
    flist_clear(&t->ins);
    for (i = 0; i < 4; i++)
        flist_clear(t->pl + i);
}

/* lower the cutoff: removal stopped paying below lc */
static void
thread_cut(hyp_thread * t, slong lc)
{
    slong i;
    if (lc >= t->lc)
        return;
    for (i = lc; i < t->nlev; i++)
    {
        pend_free(t->pq + i);
        pend_free(t->pr + i);
    }
    t->lc = lc;
}

/* after a job's subtree below level lev: its entries for the levels
   >= lev belong to the ancestors (and its cutoff applies to them) */
static void
thread_join(hyp_thread * t, hyp_thread * u, slong lev)
{
    slong i;
    thread_cut(t, u->lc);
    for (i = lev; i < t->lc; i++)
    {
        pend_append(t->pq + i, u->pq + i);
        pend_append(t->pr + i, u->pr + i);
    }
}

/* the largest prime that can cancel at a node of 2^lev leaves */
static inline double
_hyp_lev_bound(const hyp_sieve_info * I, slong lev)
{
    return I->lb[FLINT_MIN(lev, 62)];
}

/* the first level >= 1 at which entries of p can cancel, by the bounds
   lb of their factor (a row of lbq or lbr) or of any factor (lb of I);
   HYP_NEVER when none applies (no factor of the other side) */
#define HYP_NEVER 63
/* the last slot of a row of bounds holds 1 / (lb[1] - lb[0]) */
#define HYP_LB_INV 63

static inline slong
_hyp_level_row(const hyp_sieve_info * I, const double * lb, uint32_t p)
{
    slong lev;
    double x = (double) p, A, y;
    if (x <= I->cmaxp || x <= lb[1])
        return 1;
    if (!(x <= lb[HYP_NEVER - 1]))
        return HYP_NEVER;
    /* lb[t] = A 2^t + w: the least t with 2^t >= (x - w) / A, from the
       exponent of the quotient, then corrected by a step at most */
    A = lb[1] - lb[0];
    y = (x - (lb[0] - A)) * lb[HYP_LB_INV];
    {
        /* the exponent of y >= 1 (x > lb[1]), read off its bits rather
           than by ilogb, a library call */
        union { double d; uint64_t u; } c;
        c.d = y;
        lev = (slong) ((c.u >> 52) & 0x7ff) - 1023 + 1;
    }
    lev = FLINT_MAX(lev, 2);
    while (lev > 2 && x <= lb[lev - 1])
        lev--;
    while (lev < HYP_NEVER && x > lb[lev])
        lev++;
    return lev;
}

static inline slong
_hyp_level(const hyp_sieve_info * I, uint32_t p)
{
    if (I->nobound)
        return 1;
    return _hyp_level_row(I, I->lb, p);
}

static int
_fe_cmp(const void * a, const void * b)
{
    uint32_t x = ((const hyp_fe *) a)->p, y = ((const hyp_fe *) b)->p;
    return (x > y) - (x < y);
}

/* the lists of leaf j, a node at level lev (0 unless it ends a range
   short of its level): the entries of the levels <= max(lev, 1) in the
   list, the others deferred */
static void
_seg_take(hyp_flist * l, hyp_pend * pend, hyp_thread * ts,
    const hyp_sieve_info * I, const hyp_fe * src, slong n,
    const hyp_fe * cof, slong nc, slong j, slong lev)
{
    slong i, k, lv = FLINT_MAX(lev, 1);
    double b = I->nobound ? 1e300 : FLINT_MAX(_hyp_lev_bound(I, lv),
        I->cmaxp);
    hyp_flist * x = &ts->ins;

    /* the active prefix of the sorted list */
    for (k = 0; k < n && (double) src[k].p <= b; k++)
        ;
    flist_fit(l, k);
    memcpy(l->d, src, k * sizeof(hyp_fe));
    l->len = k;
    l->ok = 1;

    for (i = k; i < n; i++)
    {
        slong v = _hyp_level(I, src[i].p);
        if (v < ts->lc)
            pend_push(pend + v, (uint32_t) j, src[i].p, src[i].e);
    }

    /* the cofactors: deferred, or (a leaf above level 0) sorted into the
       list */
    x->len = 0;
    for (i = 0; i < nc; i++)
    {
        slong v = cof[i].e >> 24;
        uint32_t e = cof[i].e & 0xffffff;
        if (v <= lev)
        {
            flist_fit(x, x->len + 1);
            x->d[x->len].p = cof[i].p;
            x->d[x->len++].e = e;
        }
        else if (v < ts->lc)
            pend_push(pend + v, (uint32_t) j, cof[i].p, e);
    }
    if (x->len != 0)
    {
        slong a, w;
        qsort(x->d, x->len, sizeof(hyp_fe), _fe_cmp);
        for (a = 1, w = 1; a < x->len; a++)
        {
            if (x->d[a].p == x->d[w - 1].p)
                x->d[w - 1].e = _esat((uint64_t) x->d[w - 1].e + x->d[a].e);
            else
                x->d[w++] = x->d[a];
        }
        x->len = w;
        x->ok = 1;
        flist_add(l, x, &ts->scr);
    }
}

static void
seg_get(hyp_flist * fq, hyp_flist * fr, int need_r, hyp_thread * ts,
    const hyp_sieve_info * I, slong j, slong lev)
{
    hyp_sieve_seg * S = &ts->seg;
    slong l;

    if (j < S->j0 || j >= S->j1)
        seg_fill(S, I, j);

    l = j - S->j0;
    _seg_take(fq, ts->pq, ts, I, S->sQ + l * I->cQ, S->nQ[l],
        S->fQ + l * I->dQ, S->mQ[l], j, lev);

    if (need_r)
        _seg_take(fr, ts->pr, ts, I, S->sR + l * I->cR, S->nR[l],
            S->fR + l * I->dR, S->mR[l], j, lev);
    else
    {
        fr->ok = 0;
        fr->len = 0;
    }
}

/* the deferred entries of this level (leaves of the node, split at leaf
   jm) as sorted lists for the left and right child */
static void
_pend_split(hyp_flist * left, hyp_flist * right, hyp_pend * P, slong jm)
{
    slong i;

    left->len = right->len = 0;
    left->ok = right->ok = 1;
    if (P->len == 0)
        return;

    _ev_radix_sort_p(P->d, 0, P->len);
    flist_fit(left, P->len);
    flist_fit(right, P->len);
    for (i = 0; i < P->len; i++)
    {
        hyp_flist * t = ((slong) P->d[3 * i] >= jm) ? right : left;
        uint32_t p = P->d[3 * i + 1], e = P->d[3 * i + 2];
        if (t->len > 0 && t->d[t->len - 1].p == p)
            t->d[t->len - 1].e = _esat((uint64_t) t->d[t->len - 1].e + e);
        else
        {
            t->d[t->len].p = p;
            t->d[t->len++].e = e;
        }
    }
    P->len = 0;
}

/* the pending lists of a merge: Q left, Q right, R left, R right */
static void
hyp_pending(hyp_flist * pl, slong lev, slong jm, hyp_thread * ts)
{
    if (lev >= ts->nlev)
    {
        pl[0].len = pl[1].len = pl[2].len = pl[3].len = 0;
        pl[0].ok = pl[1].ok = pl[2].ok = pl[3].ok = 1;
        return;
    }
    _pend_split(pl + 0, pl + 1, ts->pq + lev, jm);
    _pend_split(pl + 2, pl + 3, ts->pr + lev, jm);
}

/* partial factorization of |c| into l (sorted): complete when c is a
   word, else trial division by the sieve primes */
static void
_content_factor(hyp_flist * l, const fmpz_t c, const hyp_sieve_info * I)
{
    slong i;

    l->len = 0;
    l->ok = 1;

    if (fmpz_abs_fits_ui(c))
    {
        n_factor_t fac;
        ulong a;
        fmpz_t t;
        fmpz_init(t);
        fmpz_abs(t, c);
        a = fmpz_get_ui(t);
        fmpz_clear(t);
        if (a <= 1)
            return;
        n_factor_init(&fac);
        n_factor(&fac, a, 0);
        flist_fit(l, fac.num);
        for (i = 0; i < fac.num; i++)
            if (fac.p[i] <= UINT32_MAX)
            {
                l->d[l->len].p = (uint32_t) fac.p[i];
                l->d[l->len++].e = (uint32_t) fac.exp[i];
            }
    }
    else
    {
        fmpz_t t, p;
        fmpz_init(t);
        fmpz_init(p);
        fmpz_abs(t, c);
        for (i = 0; i < I->np; i++)
        {
            slong e;
            fmpz_set_ui(p, I->primes[i]);
            e = fmpz_remove(t, t, p);
            if (e > 0)
            {
                flist_fit(l, l->len + 1);
                l->d[l->len].p = I->primes[i];
                l->d[l->len++].e = _esat(e);
            }
        }
        fmpz_clear(t);
        fmpz_clear(p);
    }
}

static void
sieve_info_clear(hyp_sieve_info * I)
{
    flint_free(I->lbq);
    flint_free(I->u);
    flint_free(I->primes);
    flist_clear(&I->cq);
    flist_clear(&I->cr);
}

/* The automatic choice.  The removal can take out at most the tracked
   bits of R per term (and of Q), and pays in proportion to their share
   r of the bits Q gains per term, against an overhead per term (the
   sieve and the lists) that only the larger precisions amortize.  The
   measured break-even: about 3-5 10^5 digits for zeta(3) (r = 0.85)
   and Catalan's constant (0.91), 3-6 10^5 for log 2 (0.75), 1-1.5 10^6
   for pi (0.52); Zuniga's log series (about 0.3) were not seen to gain
   below 2 10^6. */
static slong
_hyp_gcd_min_limbs(double r)
{
    if (r >= 0.8)
        return 24576;
    if (r >= 0.65)
        return 32768;
    if (r >= 0.45)
        return 98304;
    return WORD_MAX;
}

static double
_log2_abs_fmpz(const fmpz_t c)
{
    slong e;
    double m;
    if (fmpz_is_zero(c))
        return 0.0;
    m = fmpz_get_d_2exp(&e, c);
    return log2(fabs(m)) + (double) e;
}

/* the share r of the bits of Q per term (at the middle term of N) that
   the tracked bits of R and Q could cancel: their linear factors, and
   their contents when inc_cr resp. inc_cq */
static double
_hyp_gcd_share(const fmpz * Q, slong Qlen, const fmpz * R, slong Rlen,
    slong N, int inc_cq, int inc_cr)
{
    fmpz_poly_t f;
    fmpz_poly_factor_t fac;
    fmpz_t kk, val;
    double km = (double) FLINT_MAX(N / 2, 1), b[2] = { 0.0, 0.0 }, qt;
    slong i, j;

    if (Rlen == 0)
        return 0.0;

    fmpz_poly_init(f);
    fmpz_init_set_ui(kk, (ulong) km);
    fmpz_init(val);
    _fmpz_poly_evaluate_fmpz(val, Q, Qlen, kk);
    qt = _log2_abs_fmpz(val);

    for (j = 0; j < 2; j++)
    {
        const fmpz * F = j ? R : Q;
        slong len = j ? Rlen : Qlen;
        fmpz_poly_factor_init(fac);
        fmpz_poly_fit_length(f, len);
        _fmpz_vec_set(f->coeffs, F, len);
        _fmpz_poly_set_length(f, len);
        _fmpz_poly_normalise(f);
        fmpz_poly_factor(fac, f);
        for (i = 0; i < fac->num; i++)
            if (fac->p[i].length == 2)
            {
                fmpz_mul_ui(val, fac->p[i].coeffs + 1, (ulong) km);
                fmpz_add(val, val, fac->p[i].coeffs);
                b[j] += fac->exp[i] * _log2_abs_fmpz(val);
            }
        if (j ? inc_cr : inc_cq)
            b[j] += _log2_abs_fmpz(&fac->c);
        fmpz_poly_factor_clear(fac);
    }

    fmpz_poly_clear(f);
    fmpz_clear(kk);
    fmpz_clear(val);
    return (qt > 0.0) ? FLINT_MIN(b[0], b[1]) / qt : 0.0;
}

/* set up the sieve for Q and R (as tracked: the contents of Q resp. R
   included when inc_cq resp. inc_cr); returns 0 when there is nothing
   to gain or the values exceed the supported range */
static int
sieve_info_init(hyp_sieve_info * I, const fmpz * Q, slong Qlen,
    const fmpz * R, slong Rlen, int inc_cq, int inc_cr, slong N, slong L)
{
    fmpz_poly_t f;
    fmpz_poly_factor_t fq, fr;
    slong i, j, nf, B, cap;
    double M = 0.0;
    int ok = 0;
    fmpz_t cq, cr;

    memset(I, 0, sizeof(hyp_sieve_info));
    flist_init(&I->cq);
    flist_init(&I->cr);
    if (Rlen == 0)
        return 0;

    fmpz_poly_factor_init(fq);
    fmpz_poly_factor_init(fr);
    fmpz_poly_init(f);
    fmpz_init(cq);
    fmpz_init(cr);

    fmpz_poly_fit_length(f, Qlen);
    _fmpz_vec_set(f->coeffs, Q, Qlen);
    _fmpz_poly_set_length(f, Qlen);
    fmpz_poly_factor(fq, f);
    fmpz_poly_fit_length(f, Rlen);
    _fmpz_vec_set(f->coeffs, R, Rlen);
    _fmpz_poly_set_length(f, Rlen);
    _fmpz_poly_normalise(f);
    fmpz_poly_factor(fr, f);
    fmpz_set(cq, &fq->c);
    fmpz_set(cr, &fr->c);

    /* the distinct linear factors */
    cap = fq->num + fr->num;
    I->u = flint_malloc((cap + 1) * (2 * sizeof(ulong) + 2 * sizeof(uint32_t)));
    I->v = (slong *) (I->u + cap + 1);
    I->mQ = (uint32_t *) (I->v + cap + 1);
    I->mR = I->mQ + cap + 1;
    nf = 0;
    for (j = 0; j < 2; j++)
    {
        fmpz_poly_factor_struct * F = j ? fr : fq;
        for (i = 0; i < F->num; i++)
        {
            const fmpz * c = F->p[i].coeffs;
            slong t;
            if (F->p[i].length != 2)
                continue;
            /* primitive, positive leading coefficient */
            if (!fmpz_fits_si(c + 1) || !fmpz_fits_si(c)
                || fmpz_bits(c + 1) > FLINT_BITS - 4
                || fmpz_bits(c) > FLINT_BITS - 4)
                continue;
            for (t = 0; t < nf; t++)
                if (I->u[t] == fmpz_get_ui(c + 1)
                    && I->v[t] == fmpz_get_si(c))
                    break;
            if (t == nf)
            {
                I->u[nf] = fmpz_get_ui(c + 1);
                I->v[nf] = fmpz_get_si(c);
                I->mQ[nf] = I->mR[nf] = 0;
                nf++;
            }
            if (j)
                I->mR[t] += F->exp[i];
            else
                I->mQ[t] += F->exp[i];
        }
    }
    I->nf = nf;

    /* the largest value, and no zero values (a terminating series) */
    for (i = 0; i < nf; i++)
    {
        double m = fabs((double) I->u[i] * (double) N + (double) I->v[i]);
        M = FLINT_MAX(M, m);
        M = FLINT_MAX(M, fabs((double) I->u[i] + (double) I->v[i]));
        if (I->v[i] < 0 && (ulong) (-I->v[i]) % I->u[i] == 0
            && (ulong) (-I->v[i]) / I->u[i] <= (ulong) N)
            goto cleanup;
    }
    if (M > 0x1p60)
        goto cleanup;

    /* the matching bound and the pruning thresholds */
    {
        double qmax = 0.0, rmax = 0.0;
        I->cmul = I->cadd = 0.0;
        I->nobound = 0;
        for (i = 0; i < nf; i++)
        {
            double m = FLINT_MAX(fabs((double) I->u[i] * (double) N
                + (double) I->v[i]), fabs((double) I->u[i]
                + (double) I->v[i]));
            if (I->mQ[i])
                qmax = FLINT_MAX(qmax, m);
            if (I->mR[i])
                rmax = FLINT_MAX(rmax, m);
        }
        for (i = 0; i < nf; i++)
            for (j = 0; j < nf; j++)
                if (I->mQ[i] && I->mR[j])
                {
                    /* u_Q u_R d + (u_R v_Q - u_Q v_R), Q factor i, R
                       factor j, d = k - j >= 1 */
                    double a = (double) I->u[i] * (double) I->u[j];
                    double w = (double) I->u[j] * (double) I->v[i]
                        - (double) I->u[i] * (double) I->v[j];
                    if (w <= -a && fmod(-w, a) == 0.0)
                        I->nobound = 1;
                    I->cmul = FLINT_MAX(I->cmul, a);
                    I->cadd = FLINT_MAX(I->cadd, fabs(w));
                }
        /* entries are kept below 2^32 */
        I->qprune = (ulong) FLINT_MIN(rmax, (double) UINT32_MAX);
        I->rprune = (ulong) FLINT_MIN(qmax, (double) UINT32_MAX);
    }

    /* multiplicities are packed with a level in the cofactor entries */
    for (i = 0; i < nf; i++)
        if (I->mQ[i] >= (1u << 24) || I->mR[i] >= (1u << 24))
            goto cleanup;

    /* something to remove: a tracked factor on both sides */
    {
        int hq = 0, hr = 0;
        for (i = 0; i < nf; i++)
        {
            hq |= (I->mQ[i] != 0);
            hr |= (I->mR[i] != 0);
        }
        hq |= inc_cq && !fmpz_is_pm1(cq);
        hr |= inc_cr && !fmpz_is_pm1(cr);
        if (!hq || !hr)
            goto cleanup;
    }

    /* the primes up to sqrt(M), with their roots; at most max(4096, 16 N)
       and 2^24 (factors with huge coefficients): a cofactor may then be
       composite, which only makes it an opaque label (entries are
       disjoint parts of the values, so equal labels still divide both
       sides, and the bounds on shared divisors hold for any of them) */
    B = (slong) n_sqrt((ulong) M) + 1;
    B = FLINT_MIN(B, FLINT_MAX(4096, 16 * N));
    B = FLINT_MIN(B, WORD(1) << 24);
    {
        char * comp = flint_calloc(B + 2, 1);
        slong np = 0, pw;
        for (i = 2; i <= B; i++)
            if (!comp[i])
            {
                np++;
                for (j = i * i; j <= B; j += i)
                    comp[j] = 1;
            }
        /* primes (uint32), then pinv, plim (words), then the roots */
        pw = (np * sizeof(uint32_t) + sizeof(ulong) - 1)
            / sizeof(ulong);
        I->np = np;
        I->primes = flint_malloc((pw + 2 * np) * sizeof(ulong)
            + np * nf * sizeof(uint32_t) + 1);
        I->pinv = (ulong *) I->primes + pw;
        I->plim = I->pinv + np;
        I->roots = (uint32_t *) (I->plim + np);
        np = 0;
        for (i = 2; i <= B; i++)
            if (!comp[i])
                I->primes[np++] = (uint32_t) i;
        flint_free(comp);
    }
    for (i = 0; i < I->np; i++)
    {
        ulong p = I->primes[i], inv = 1;
        /* p^-1 mod 2^FLINT_BITS by Newton (p odd) */
        if (p != 2)
        {
            inv = p;
            for (j = 0; j < 6; j++)
                inv *= 2 - p * inv;
        }
        I->pinv[i] = inv;
        I->plim[i] = UWORD_MAX / p;
        for (j = 0; j < nf; j++)
        {
            ulong um = I->u[j] % p, vm;
            if (um == 0)
            {
                I->roots[i * nf + j] = HYP_NO_ROOT;
                continue;
            }
            vm = (I->v[j] >= 0) ? (ulong) I->v[j] % p
                : (p - (ulong) (-I->v[j]) % p) % p;
            /* root of u k + v: k = -v / u mod p */
            I->roots[i * nf + j] = (uint32_t) n_mulmod2((p - vm) % p,
                n_invmod(um, p), p);
        }
    }

    if (inc_cq)
        _content_factor(&I->cq, cq, I);
    else
        I->cq.ok = 1;
    if (inc_cr)
        _content_factor(&I->cr, cr, I);
    else
        I->cr.ok = 1;

    /* content primes cancel at any distance */
    I->cmaxp = 0.0;
    if (I->cq.len)
    {
        I->cmaxp = FLINT_MAX(I->cmaxp, (double) I->cq.d[I->cq.len - 1].p);
        I->rprune = FLINT_MAX(I->rprune, I->cq.d[I->cq.len - 1].p);
    }
    if (I->cr.len)
    {
        I->cmaxp = FLINT_MAX(I->cmaxp, (double) I->cr.d[I->cr.len - 1].p);
        I->qprune = FLINT_MAX(I->qprune, I->cr.d[I->cr.len - 1].p);
    }

    I->L = L;
    I->N = N;
    I->nleaves = (N + L - 1) / L;
    /* Per leaf: the distinct sieve primes, about E = sum over the primes
       p of 1 - (1 - min(1, L/p))^f (f the factors with a root mod p),
       given 1.5 E + 64 slots (at most every sieve prime), the content
       primes, and one cofactor per value.  Slots sized for every sieve
       prime would make a window several times larger than the lists
       it holds (and shorter, for a given memory). */
    {
        slong nq = 0, nr = 0, per, want, t;
        double eq = 0.0, er = 0.0;
        for (i = 0; i < nf; i++)
        {
            nq += (I->mQ[i] != 0);
            nr += (I->mR[i] != 0);
        }
        for (t = 0; t < I->np; t++)
        {
            double x = FLINT_MIN(1.0, (double) L / I->primes[t]);
            slong fq = 0, fr = 0;
            for (i = 0; i < nf; i++)
                if (I->roots[t * nf + i] != HYP_NO_ROOT)
                {
                    fq += (I->mQ[i] != 0);
                    fr += (I->mR[i] != 0);
                }
            eq += 1.0 - pow(1.0 - x, (double) fq);
            er += 1.0 - pow(1.0 - x, (double) fr);
        }
        I->cQ = FLINT_MIN(I->np, (slong) (1.5 * eq) + 64) + I->cq.len + 1;
        I->cR = FLINT_MIN(I->np, (slong) (1.5 * er) + 64) + I->cr.len + 1;
        I->dQ = L * nq + 1;
        I->dR = L * nr + 1;
        /* windows long enough that the roots (np nf per window) cost a
           few operations per term, within about 4 MB (1 MB was measured
           10-20% slower at 10^8 digits, 8 MB no faster) */
        per = (I->cQ + I->cR + I->dQ + I->dR) * sizeof(hyp_fe)
            + nf * L * sizeof(ulong) + 4 * sizeof(uint32_t);
        want = FLINT_MAX(64 * L, I->np * nf / 4);
        I->segleaves = (want + L - 1) / L;
        I->segleaves = FLINT_MIN(I->segleaves, (WORD(1) << 22) / per);
        I->segleaves = FLINT_MAX(I->segleaves, 8);
    }
    for (i = 0; i < 63; i++)
        I->lb[i] = I->cmul * ldexp((double) L, (int) i) + I->cadd;
    I->lb[63] = 1.0 / FLINT_MAX(I->lb[1] - I->lb[0], 1e-300);

    /* per factor: against the factors of the other side */
    I->lbq = flint_malloc(2 * 64 * (nf + 1) * sizeof(double));
    I->lbr = I->lbq + 64 * (nf + 1);
    for (i = 0; i < nf; i++)
    {
        double aq = 0.0, wq = 0.0, ar = 0.0, wr = 0.0;
        int nbq = 0, nbr = 0;
        slong t;
        for (j = 0; j < nf; j++)
        {
            double a = (double) I->u[i] * (double) I->u[j];
            if (I->mQ[i] && I->mR[j])
            {
                double w = (double) I->u[j] * (double) I->v[i]
                    - (double) I->u[i] * (double) I->v[j];
                nbq |= (w <= -a && fmod(-w, a) == 0.0);
                aq = FLINT_MAX(aq, a);
                wq = FLINT_MAX(wq, fabs(w));
            }
            if (I->mR[i] && I->mQ[j])
            {
                double w = (double) I->u[i] * (double) I->v[j]
                    - (double) I->u[j] * (double) I->v[i];
                nbr |= (w <= -a && fmod(-w, a) == 0.0);
                ar = FLINT_MAX(ar, a);
                wr = FLINT_MAX(wr, fabs(w));
            }
        }
        for (t = 0; t < 63; t++)
        {
            I->lbq[64 * i + t] = nbq ? 1e300
                : aq * ldexp((double) L, (int) t) + wq;
            I->lbr[64 * i + t] = nbr ? 1e300
                : ar * ldexp((double) L, (int) t) + wr;
        }
        I->lbq[64 * i + 63] = 1.0 / FLINT_MAX(aq * L, 1e-300);
        I->lbr[64 * i + 63] = 1.0 / FLINT_MAX(ar * L, 1e-300);
    }
    ok = 1;

cleanup:
    if (!ok)
    {
        sieve_info_clear(I);
        memset(I, 0, sizeof(hyp_sieve_info));
    }
    fmpz_poly_factor_clear(fq);
    fmpz_poly_factor_clear(fr);
    fmpz_poly_clear(f);
    fmpz_clear(cq);
    fmpz_clear(cr);
    return ok;
}

/* ---- the series ---- */

typedef struct
{
    hpoly_t P, Q, R;        /* the polynomials */
    hpoly_t Qc, Rc;         /* Q / cQ, R / cR (content mode) */
    hpoly_t A;              /* P / R (factored mode) */
    int factored;           /* R divides P */
    fmpz_t cQ, cR;          /* contents (positive) */
    int content;            /* content mode */
    int have_cR;            /* |cR| > 1 in content mode */
    slong L, J, wp;
    nn_ptr buf;
    slong bufn;
    slong vlimbs;           /* room for the polynomial values */
    mp_real_struct * cQpow;   /* cQ^(L 2^i), i < J */
    mp_real_struct * cRpow;
    mp_real_struct * tmp;     /* 5 per level */
    double tbits;           /* estimated bits per term of T, Q */
    int gcd;                /* content removal */
    const hyp_sieve_info * si;
    hyp_thread * ts;        /* this thread's sieve window and buffers */
    hyp_flist * ftmp;       /* 3 per level: fq2, fr2, g */
}
hyp_struct;

#define NBUF 11

/* the leaf when every value is a single word (and no content powers):
   the common case, with the per-term work reduced to the mpn passes */
static void
leaf_w1(mp_real_t T, mp_real_t Qo, mp_real_t Ro, slong a, slong b, int need_r,
    int content, const hpoly_struct * Rp, const hyp_struct * H)
{
    slong bufn = H->bufn, m = b - a, i;
    nn_ptr Tb = H->buf, Xb = Tb + bufn, Qs = Xb + bufn, Rs = Qs + bufn,
        Cs = Rs + bufn;
    slong tn = 0, qn = 1, rn = 1, cn = 1;
    int tneg = 0, qneg = 0, rneg = 0, cneg = 0, factored = H->factored;
    ulong xb[HYP_LEAF_TERMS_MAX], qb[HYP_LEAF_TERMS_MAX],
        rb[HYP_LEAF_TERMS_MAX], cb[HYP_LEAF_TERMS_MAX],
        rpb[HYP_LEAF_TERMS_MAX], cy;
    const ulong * rob = rb;

    FLINT_ASSERT(m <= HYP_LEAF_TERMS_MAX);

    hpoly_eval_block(xb, factored ? H->A : H->P, a, m);
    hpoly_eval_block(qb, H->Q, a, m);
    hpoly_eval_block(rb, H->R, a, m);
    if (content)
        hpoly_eval_block(cb, H->Qc, a, m);
    if (need_r && Rp != H->R)
    {
        hpoly_eval_block(rpb, Rp, a, m);
        rob = rpb;
    }

    Qs[0] = 1;
    Rs[0] = 1;
    Cs[0] = 1;

    for (i = m - 1; i >= 0; i--)
    {
        slong xv = (slong) xb[i], qv = (slong) qb[i], rv = (slong) rb[i];
        ulong x = FLINT_UABS(xv), q = FLINT_UABS(qv), r = FLINT_UABS(rv);
        int xneg = (xv < 0) ^ qneg, rs = (rv < 0);

        if (factored)
        {
            /* T <- R(k) (T + A(k) Qs) */
            if (x != 0)
                _addmul_1_signed(Tb, &tn, &tneg, Qs, qn, x, xneg);
            if (r == 0)
            {
                tn = 0;
            }
            else if (tn != 0)
            {
                cy = mpn_mul_1(Tb, Tb, tn, r);
                Tb[tn] = cy;
                tn += (cy != 0);
                tneg ^= rs;
            }
        }
        else
        {
            /* T <- P(k) Qs + R(k) T */
            slong xn;

            if (x == 1)
            {
                flint_mpn_copyi(Xb, Qs, qn);
                xn = qn;
            }
            else if (x == 0)
            {
                xn = 0;
            }
            else
            {
                cy = mpn_mul_1(Xb, Qs, qn, x);
                Xb[qn] = cy;
                xn = qn + (cy != 0);
            }
            if (tn != 0 && r != 0)
                _addmul_1_signed(Xb, &xn, &xneg, Tb, tn, r, rs ^ tneg);
            FLINT_SWAP(nn_ptr, Tb, Xb);
            tn = xn;
            tneg = xneg;
        }

        FLINT_ASSERT(q != 0);
        cy = mpn_mul_1(Qs, Qs, qn, q);
        Qs[qn] = cy;
        qn += (cy != 0);
        qneg ^= (qv < 0);

        if (content)
        {
            slong cv = (slong) cb[i];
            cy = mpn_mul_1(Cs, Cs, cn, FLINT_UABS(cv));
            Cs[cn] = cy;
            cn += (cy != 0);
            cneg ^= (cv < 0);
        }

        if (need_r)
        {
            slong ov = (slong) rob[i];
            ulong o = FLINT_UABS(ov);

            if (o == 0)
            {
                rn = 0;
            }
            else if (rn != 0)
            {
                cy = mpn_mul_1(Rs, Rs, rn, o);
                Rs[rn] = cy;
                rn += (cy != 0);
            }
            rneg ^= (ov < 0);
        }

        FLINT_ASSERT(tn + 4 <= bufn && qn + 4 <= bufn && rn + 4 <= bufn);
    }

    _mp_real_set_mpn_2exp(T, Tb, tn, 0);
    if (tneg)
        mp_real_neg(T, T);
    if (content)
    {
        _mp_real_set_mpn_2exp(Qo, Cs, cn, 0);
        if (cneg)
            mp_real_neg(Qo, Qo);
    }
    else
    {
        _mp_real_set_mpn_2exp(Qo, Qs, qn, 0);
        if (qneg)
            mp_real_neg(Qo, Qo);
    }
    if (need_r)
    {
        _mp_real_set_mpn_2exp(Ro, Rs, rn, 0);
        if (rneg)
            mp_real_neg(Ro, Ro);
    }
}

/* exact T, Q (or Q'), R (or R') over [a, b) as mp_real_t balls; single-word
   polynomial values are computed for the whole block up front */
static void
leaf(mp_real_t T, mp_real_t Qo, mp_real_t Ro, slong a, slong b, int need_r,
    int full, const hyp_struct * H)
{
    slong bufn = H->bufn, k, m = b - a;
    nn_ptr Tb = H->buf, Xb = Tb + bufn, Yb = Xb + bufn,
        Qa = Yb + bufn, Qb = Qa + bufn, Pa = Qb + bufn, Pb = Pa + bufn,
        Ra = Pb + bufn, Rb = Ra + bufn, vbuf = Rb + bufn;
    slong tn = 0, qn = 1, pn = 1, rn = 1;
    int tneg = 0, qneg = 0, pneg = 0, rneg = 0;
    int content = H->content && !full;
    int factored = H->factored;
    const hpoly_struct * X = factored ? H->A : H->P;   /* A or P */
    const hpoly_struct * Rp = (content && H->have_cR) ? H->Rc : H->R;
    ulong xb[HYP_LEAF_TERMS_MAX], qb[HYP_LEAF_TERMS_MAX],
        rb[HYP_LEAF_TERMS_MAX], qcb[HYP_LEAF_TERMS_MAX],
        rcb[HYP_LEAF_TERMS_MAX];
    nn_ptr xv = vbuf, qv = xv + FLINT_MAX(H->P->W, H->A->W),
        rv = qv + H->Q->W, qcv = rv + H->R->W, rcv = qcv + H->Q->W;

    FLINT_ASSERT(m <= HYP_LEAF_TERMS_MAX);

    if (X->W == 1 && H->Q->W == 1 && H->R->W == 1
        && (!content || H->Qc->W == 1) && (!need_r || Rp->W == 1))
    {
        leaf_w1(T, Qo, Ro, a, b, need_r, content, Rp, H);
        return;
    }

    /* block values, indexed by k - a */
    if (X->W == 1)
        hpoly_eval_block(xb, X, a, m);
    if (H->Q->W == 1)
        hpoly_eval_block(qb, H->Q, a, m);
    if (H->R->W == 1)
        hpoly_eval_block(rb, H->R, a, m);
    if (content && H->Qc->W == 1)
        hpoly_eval_block(qcb, H->Qc, a, m);
    if (need_r && Rp != H->R && Rp->W == 1)
        hpoly_eval_block(rcb, Rp, a, m);

    Qa[0] = 1;      /* Qs: the full Q product (T needs it) */
    Pa[0] = 1;      /* Q' product (content mode) */
    Ra[0] = 1;      /* R or R' product */

    for (k = b - 1; k >= a; k--)
    {
        slong i = k - a, xvn, qvn, rvn;
        int xs, qs, rs;

        xs = hpoly_value(xv, &xvn, X, xb, i, k);
        qs = hpoly_value(qv, &qvn, H->Q, qb, i, k);
        rs = hpoly_value(rv, &rvn, H->R, rb, i, k);

        if (factored)
        {
            /* T <- R(k) (T + A(k) Qs) */
            if (xvn == 1)
            {
                _addmul_1_signed(Tb, &tn, &tneg, Qa, qn, xv[0], xs ^ qneg);
            }
            else if (xvn != 0)
            {
                slong yn = _mul_small(Yb, Qa, qn, xv, xvn);
                _add_signed(Tb, &tn, &tneg, Yb, yn, xs ^ qneg);
            }

            if (rvn == 0)
            {
                tn = 0;
            }
            else if (tn != 0)
            {
                if (rvn == 1)
                {
                    tn = _mul_small(Tb, Tb, tn, rv, 1);
                }
                else
                {
                    tn = _mul_small(Xb, Tb, tn, rv, rvn);
                    FLINT_SWAP(nn_ptr, Tb, Xb);
                }
                tneg ^= rs;
            }
        }
        else
        {
            /* T <- P(k) Qs + R(k) T */
            slong xn;
            int xneg = xs ^ qneg;

            if (xvn == 1 && xv[0] == 1)
            {
                flint_mpn_copyi(Xb, Qa, qn);
                xn = qn;
            }
            else
                xn = (xvn != 0) ? _mul_small(Xb, Qa, qn, xv, xvn) : 0;

            if (tn != 0 && rvn != 0)
            {
                if (rvn == 1)
                {
                    _addmul_1_signed(Xb, &xn, &xneg, Tb, tn, rv[0],
                        rs ^ tneg);
                }
                else
                {
                    slong yn = _mul_small(Yb, Tb, tn, rv, rvn);
                    _add_signed(Xb, &xn, &xneg, Yb, yn, rs ^ tneg);
                }
            }
            FLINT_SWAP(nn_ptr, Tb, Xb);
            tn = xn;
            tneg = xneg;
        }

        /* Qs <- Q(k) Qs */
        FLINT_ASSERT(qvn != 0);
        if (qvn == 1)
        {
            qn = _mul_small(Qa, Qa, qn, qv, 1);
        }
        else
        {
            qn = _mul_small(Qb, Qa, qn, qv, qvn);
            FLINT_SWAP(nn_ptr, Qa, Qb);
        }
        qneg ^= qs;

        if (content)
        {
            slong cvn;
            int cs = hpoly_value(qcv, &cvn, H->Qc, qcb, i, k);
            if (cvn == 1)
            {
                pn = _mul_small(Pa, Pa, pn, qcv, 1);
            }
            else
            {
                pn = _mul_small(Pb, Pa, pn, qcv, cvn);
                FLINT_SWAP(nn_ptr, Pa, Pb);
            }
            pneg ^= cs;
        }

        if (need_r)
        {
            slong cvn;
            int cs;
            nn_ptr cv;

            if (Rp == H->R)
            {
                cv = rv;
                cvn = rvn;
                cs = rs;
            }
            else
            {
                cv = rcv;
                cs = hpoly_value(rcv, &cvn, Rp, rcb, i, k);
            }

            if (cvn == 0)
                rn = 0;
            else if (cvn == 1)
                rn = _mul_small(Ra, Ra, rn, cv, 1);
            else if (rn != 0)
            {
                rn = _mul_small(Rb, Ra, rn, cv, cvn);
                FLINT_SWAP(nn_ptr, Ra, Rb);
            }
            rneg ^= cs;
        }

        FLINT_ASSERT(tn + 4 <= bufn && qn + 4 <= bufn && rn + 4 <= bufn
            && pn + 4 <= bufn);
    }

    _mp_real_set_mpn_2exp(T, Tb, tn, 0);
    if (tneg)
        mp_real_neg(T, T);

    if (content)
    {
        _mp_real_set_mpn_2exp(Qo, Pa, pn, 0);
        if (pneg)
            mp_real_neg(Qo, Qo);
    }
    else
    {
        _mp_real_set_mpn_2exp(Qo, Qa, qn, 0);
        if (qneg)
            mp_real_neg(Qo, Qo);
    }

    if (need_r)
    {
        _mp_real_set_mpn_2exp(Ro, Ra, rn, 0);
        if (rneg)
            mp_real_neg(Ro, Ro);
    }
}

/* T = T1 Q2 + R1 T2 (with the content powers pq, pr when non-NULL),
   Q = Q1 Q2, R = R1 R2; (T, Q, R) hold the left values on entry;
   (T2, Q2, R2, U, V) is scratch */
static void
merge(mp_real_t T, mp_real_t Q, mp_real_t R, mp_real_t T2, mp_real_t Q2, mp_real_t R2,
    mp_real_t U, mp_real_t V, const mp_real_struct * pq, const mp_real_struct * pr,
    int need_r, slong wp)
{
    mp_real_mul(U, T, Q2, wp);
    if (pq != NULL)
    {
        mp_real_mul(V, U, pq, wp);
        mp_real_swap(U, V);
    }

    mp_real_mul(V, R, T2, wp);
    if (pr != NULL)
    {
        mp_real_mul(T2, V, pr, wp);
        mp_real_swap(T2, V);
    }

    mp_real_add(T, U, V, wp);

    mp_real_mul(U, Q, Q2, wp);
    mp_real_swap(Q, U);

    if (need_r)
    {
        mp_real_mul(U, R, R2, wp);
        mp_real_swap(R, U);
    }
}

/* the content removal at a node: g = gcd(R1, Q2) read off the lists,
   R1 and Q2 divided by it (T = T1 Q2 + R1 T2, Q = Q1 Q2 and R = R1 R2
   all have the factor g, and the node only represents the ratios T/Q,
   R/Q), when both are exact integers */
typedef struct
{
    mp_real_struct * x;
    nn_srcptr G;
    slong gn;
    const flint_mpn_divexact_preinv_struct * pre;
}
hyp_div_struct;

static void
_hyp_div_job(void * arg)
{
    hyp_div_struct * D = (hyp_div_struct *) arg;
    _mp_real_divexact_mpn_pre(D->x, D->G, D->gn, D->pre);
}

/* from this size of G, its 2-adic inverse is computed once for both
   divisions (measured 10-30% faster beyond a few thousand limbs, no
   gain below) */
#define HYP_PREINV_LIMBS 1000

/* The removal pays when the bits it takes out of the products of this
   node and its h ancestors (about four products each, with one operand
   g bits shorter) exceed the cost of the two exact divisions.  A 2-adic
   division costs in proportion to its quotient (n - g limbs for n by g;
   near the root R1 is often nearly all g, and its division nearly
   free), plus an inverse of about the divisor's size, and per limb
   about 0.15-0.4 products of the same size below a few hundred limbs
   (a short divisor makes it nearly linear) and 1.5-2.4 beyond a few
   thousand (measured); rho is that, doubled (a product of two n-limb
   operands counted as 2n), times 1.5, the best of the factors tried
   (1, 1.5, 2, 2.5) at 10^5 to 3 10^6 digits.  Counting the quotients
   rather than the dividends was measured 2-3% faster at 10^7 digits.
   Returns 0 when the removal does not pay here, and the lists of the
   node can then be dropped: the ancestors have larger operands and
   fewer levels above. */
static int
_hyp_remove_pays(double gbits, slong qn, slong rn, slong h)
{
    double gl = gbits / FLINT_BITS;
    double q = FLINT_MAX((double) qn - gl, 0.0);
    double r = FLINT_MAX((double) rn - gl, 0.0);
    double n = FLINT_MAX(q, r), rho;

    if (n < 64)
        rho = 0.45;
    else if (n < 512)
        rho = 0.9;
    else if (n < 2048)
        rho = 2.25;
    else
        rho = 4.5;

    return gbits * 4.0 * (double) (h + 1)
        >= rho * (double) FLINT_BITS * (q + r + gl);
}

static int
hyp_remove(mp_real_t R1, mp_real_t Q2, hyp_flist * fr1, hyp_flist * fq2,
    hyp_flist * pl, hyp_flist * g, slong h, slong terms,
    const hyp_sieve_info * I, int par)
{
    nn_ptr G;
    slong gn;
    ulong sbuf[HYP_SMALL_G];
    double gb;

    if (!fr1->ok || !fq2->ok || !_mp_real_is_exact_int(R1)
        || !_mp_real_is_exact_int(Q2))
        return 0;
    if (fr1->len + pl[2].len == 0 || fq2->len + pl[1].len == 0)
        return 1;

    /* only primes up to the bound for this node can be shared */
    {
        double b = I->nobound ? (double) UINT32_MAX
            : FLINT_MAX(I->cmul * (double) terms + I->cadd, I->cmaxp);
        flist_gcd_remove(g, fr1, pl + 2, fq2, pl + 1,
            (ulong) FLINT_MIN(b, (double) UINT32_MAX));
    }
    if (g->len == 0)
        return 1;

    gb = flist_bits(g);
    if (!_hyp_remove_pays(gb, Q2->size, R1->size, h))
        return 0;

    gn = flist_expand(&G, g, sbuf);
    {
        flint_mpn_divexact_preinv_t pre;
        const flint_mpn_divexact_preinv_struct * pp = NULL;
        if (gn >= HYP_PREINV_LIMBS)
        {
            flint_mpn_divexact_preinv_init(pre, G, gn);
            pp = pre;
        }
        if (par)
        {
            hyp_div_struct A, B;
            A.x = R1; A.G = G; A.gn = gn; A.pre = pp;
            B.x = Q2; B.G = G; B.gn = gn; B.pre = pp;
            _mp_real_parallel_pair(_hyp_div_job, &A, _hyp_div_job, &B);
        }
        else
        {
            _mp_real_divexact_mpn_pre(R1, G, gn, pp);
            _mp_real_divexact_mpn_pre(Q2, G, gn, pp);
        }
        if (pp != NULL)
            flint_mpn_divexact_preinv_clear(pre);
    }
    if (G != sbuf)
        flint_free(G);
    return 1;
}

static void
_flist_drop(hyp_flist * l)
{
    flist_clear(l);
}

/* the lists of the merged node (Q = Q1 Q2, R = R1 R2), dropped when the
   values have become inexact; the right child's lists are emptied */
static void
hyp_lists_merge(hyp_flist * FQ, hyp_flist * FR, hyp_flist * fq2,
    hyp_flist * fr2, hyp_flist * pl, const mp_real_t Q, const mp_real_t R,
    int need_r, hyp_flist * scr)
{
    /* the pending lists are short: pl[0] + pl[1] (resp. pl[2] + pl[3])
       first, then one pass over the long lists */
    if (_mp_real_is_exact_int(Q) && FQ->ok && fq2->ok)
    {
        flist_add(pl + 0, pl + 1, scr);
        flist_add3(FQ, fq2, pl + 0, scr);
    }
    else
        _flist_drop(FQ);

    if (need_r && _mp_real_is_exact_int(R) && FR->ok && fr2->ok)
    {
        flist_add(pl + 2, pl + 3, scr);
        flist_add3(FR, fr2, pl + 2, scr);
    }
    else
        _flist_drop(FR);

    fq2->len = 0;
    fq2->ok = 0;
    fr2->len = 0;
    fr2->ok = 0;
}

/* the removal stops at this node (it does not pay, or the values are no
   longer exact): the lists are dropped here and above, and the entries
   deferred to the levels above are discarded */
static void
hyp_stop(hyp_flist * FQ, hyp_flist * FR, hyp_flist * f2, slong lev,
    hyp_thread * ts)
{
    _flist_drop(FQ);
    _flist_drop(FR);
    _flist_drop(f2);
    _flist_drop(f2 + 1);
    thread_cut(ts, lev + 1);
}

static void bsplit_pow2(mp_real_t T, mp_real_t Q, mp_real_t R, hyp_flist * FQ,
    hyp_flist * FR, slong a0, slong k, slong lev, int need_r,
    const hyp_struct * H);
static void bsplit_range(mp_real_t T, mp_real_t Q, mp_real_t R, hyp_flist * FQ,
    hyp_flist * FR, slong a, slong b, slong lev, int need_r,
    const hyp_struct * H);

/* a subtree as a job for _mp_real_parallel_pair; with own set it runs on
   private leaf buffers and level temporaries (the shared H only
   supplies the read-only polynomials and content powers), freed at the
   end */
typedef struct
{
    mp_real_struct * T, * Q, * R;
    hyp_flist * FQ, * FR;
    slong a, b, lev;
    int need_r, pow2, own;
    const hyp_struct * H;
    hyp_thread * ts;            /* own: this job's thread state */
}
hyp_job;

static void
_hyp_job(void * arg)
{
    hyp_job * J = (hyp_job *) arg;
    hyp_struct H2;
    const hyp_struct * H = J->H;
    slong i, nt = 5 * (J->lev + 1), nf = 3 * (J->lev + 1);

    if (J->own)
    {
        H2 = *J->H;
        H2.buf = flint_malloc((NBUF * H2.bufn + H2.vlimbs) * sizeof(ulong));
        H2.tmp = flint_malloc(nt * sizeof(mp_real_struct));
        for (i = 0; i < nt; i++)
            mp_real_init(H2.tmp + i);
        if (H2.gcd)
        {
            H2.ts = J->ts;
            H2.ftmp = flint_malloc(nf * sizeof(hyp_flist));
            for (i = 0; i < nf; i++)
                flist_init(H2.ftmp + i);
        }
        H = &H2;
    }

    if (J->pow2)
        bsplit_pow2(J->T, J->Q, J->R, J->FQ, J->FR, J->a, J->b, J->lev,
            J->need_r, H);
    else
        bsplit_range(J->T, J->Q, J->R, J->FQ, J->FR, J->a, J->b, J->lev,
            J->need_r, H);

    if (J->own)
    {
        for (i = 0; i < nt; i++)
            mp_real_clear(H2.tmp + i);
        flint_free(H2.tmp);
        flint_free(H2.buf);
        if (H2.gcd)
        {
            for (i = 0; i < nf; i++)
                flist_clear(H2.ftmp + i);
            flint_free(H2.ftmp);
        }
    }
}

/* subtrees of at least HYP_PAR_LEAVES leaves may run on two threads,
   by the rule of MP_REAL_PAR_CAP: 0 serial, 1 the halves on two
   threads, 2 above the cap (the merge on two threads) */
#define HYP_PAR_LEAVES 64

static int
_hyp_par(slong leaves, slong terms, const hyp_struct * H)
{
    if (leaves < HYP_PAR_LEAVES || flint_get_num_threads() < 2)
        return 0;
    return ((double) terms * H->tbits
        > MP_REAL_PAR_CAP * FLINT_BITS * (double) H->wp) ? 2 : 1;
}

static void
_hyp_halves(mp_real_struct * T, mp_real_struct * Q, mp_real_struct * R_,
    hyp_flist * FQ, hyp_flist * FR, mp_real_struct * T2, hyp_flist * f2,
    slong a, slong m, slong b,
    slong lev, int need_r, int pow2, int fork, const hyp_struct * H)
{
    hyp_flist * fq2 = f2, * fr2 = (f2 == NULL) ? NULL : f2 + 1;

    if (fork)
    {
        hyp_job L, R;
        L.T = T; L.Q = Q; L.R = R_; L.FQ = FQ; L.FR = FR; L.a = a; L.b = m;
        L.lev = lev - 1; L.need_r = 1; L.pow2 = pow2; L.own = 0; L.H = H;
        R.T = T2; R.Q = T2 + 1; R.R = T2 + 2; R.FQ = fq2; R.FR = fr2;
        R.a = pow2 ? a + m : m;
        R.b = pow2 ? m : b; R.lev = lev - 1; R.need_r = need_r;
        R.pow2 = pow2; R.own = 1; R.H = H;
        L.ts = NULL;
        R.ts = NULL;
        if (H->gcd)
        {
            R.ts = flint_malloc(sizeof(hyp_thread));
            thread_init(R.ts, H->ts->nlev);
            R.ts->lc = H->ts->lc;
        }
        _mp_real_parallel_pair(_hyp_job, &L, _hyp_job, &R);
        if (H->gcd)
        {
            thread_join(H->ts, R.ts, lev);
            thread_clear(R.ts);
            flint_free(R.ts);
        }
    }
    else if (pow2)
    {
        bsplit_pow2(T, Q, R_, FQ, FR, a, m, lev - 1, 1, H);
        bsplit_pow2(T2, T2 + 1, T2 + 2, fq2, fr2, a + m, m, lev - 1, need_r,
            H);
    }
    else
    {
        bsplit_range(T, Q, R_, FQ, FR, a, m, lev - 1, 1, H);
        bsplit_range(T2, T2 + 1, T2 + 2, fq2, fr2, m, b, lev - 1, need_r, H);
    }
}

/* merge on two threads: {U = T1 Q2 (pq), Q = Q1 Q2} and {V = R1 T2 (pr),
   R = R1 R2} */
typedef struct
{
    mp_real_struct * T, * Q, * R, * T2, * Q2, * R2, * U, * V;
    const mp_real_struct * pq, * pr;
    int need_r;
    slong wp;
}
hyp_merge_struct;

static void
_hyp_merge_x(void * arg)
{
    hyp_merge_struct * M = (hyp_merge_struct *) arg;
    mp_real_mul(M->U, M->T, M->Q2, M->wp);
    if (M->pq != NULL)
        mp_real_mul(M->U, M->U, M->pq, M->wp);
    mp_real_mul(M->Q, M->Q, M->Q2, M->wp);
}

static void
_hyp_merge_y(void * arg)
{
    hyp_merge_struct * M = (hyp_merge_struct *) arg;
    mp_real_mul(M->V, M->R, M->T2, M->wp);
    if (M->pr != NULL)
        mp_real_mul(M->V, M->V, M->pr, M->wp);
    if (M->need_r)
        mp_real_mul(M->R, M->R, M->R2, M->wp);
}

static void
merge_par(mp_real_t T, mp_real_t Q, mp_real_t R, mp_real_t T2, mp_real_t Q2,
    mp_real_t R2, mp_real_t U, mp_real_t V, const mp_real_struct * pq,
    const mp_real_struct * pr, int need_r, slong wp)
{
    hyp_merge_struct M;
    M.T = T; M.Q = Q; M.R = R; M.T2 = T2; M.Q2 = Q2; M.R2 = R2;
    M.U = U; M.V = V; M.pq = pq; M.pr = pr; M.need_r = need_r; M.wp = wp;
    _mp_real_parallel_pair(_hyp_merge_x, &M, _hyp_merge_y, &M);
    mp_real_add(T, U, V, wp);
}

/* free this level's temporaries when large: they are otherwise held
   through all the merges above */
#define HYP_KEEP 2048

static void
_hyp_release(mp_real_struct * T2, hyp_flist * f2)
{
    slong i;
    if (T2[0].alloc > HYP_KEEP || T2[1].alloc > HYP_KEEP)
        for (i = 0; i < 5; i++)
        {
            mp_real_clear(T2 + i);
            mp_real_init(T2 + i);
        }
    if (f2 != NULL)
        for (i = 0; i < 3; i++)
            if (f2[i].alloc > HYP_KEEP)
                flist_clear(f2 + i);
}

/* a leaf, with its lists */
static void
hyp_leaf(mp_real_t T, mp_real_t Q, mp_real_t R, hyp_flist * FQ,
    hyp_flist * FR, slong a, slong b, slong lev, int need_r, int full,
    const hyp_struct * H)
{
    leaf(T, Q, R, a, b, need_r, full, H);
    if (H->gcd)
        seg_get(FQ, FR, need_r, H->ts, H->si, (a - 1) / H->L, lev);
}

/* content mode: k leaves (a power of two) from leaf a0, at level
   lev = log2(k) */
static void
bsplit_pow2(mp_real_t T, mp_real_t Q, mp_real_t R, hyp_flist * FQ,
    hyp_flist * FR, slong a0, slong k, slong lev, int need_r,
    const hyp_struct * H)
{
    if (k == 1)
    {
        hyp_leaf(T, Q, R, FQ, FR, 1 + a0 * H->L, 1 + (a0 + 1) * H->L, 0,
            need_r, 0, H);
    }
    else
    {
        mp_real_struct * T2 = H->tmp + 5 * lev;
        hyp_flist * f2 = H->gcd ? H->ftmp + 3 * lev : NULL;
        slong h = k / 2;

        int par = _hyp_par(k, k * H->L, H);

        _hyp_halves(T, Q, R, FQ, FR, T2, f2, a0, h, 0, lev, need_r, 1,
            par == 1, H);
        if (H->gcd)
        {
            hyp_pending(H->ts->pl, lev, a0 + h, H->ts);
            if (!hyp_remove(R, T2 + 1, FR, f2, H->ts->pl, f2 + 2,
                    H->J - lev, k * H->L, H->si, par == 2))
                hyp_stop(FQ, FR, f2, lev, H->ts);
        }
        (par == 2 ? merge_par : merge)(T, Q, R, T2, T2 + 1, T2 + 2,
            T2 + 3, T2 + 4, H->cQpow + (lev - 1),
            H->have_cR ? H->cRpow + (lev - 1) : NULL, need_r, H->wp);
        if (H->gcd)
            hyp_lists_merge(FQ, FR, f2, f2 + 1, H->ts->pl, Q, R, need_r,
                &H->ts->scr);
        _hyp_release(T2, f2);
    }
}

/* carry mode: the terms [a, b), in leaves of L terms (the last one
   possibly shorter), at most 2^lev leaves */
static void
bsplit_range(mp_real_t T, mp_real_t Q, mp_real_t R, hyp_flist * FQ,
    hyp_flist * FR, slong a, slong b, slong lev, int need_r,
    const hyp_struct * H)
{
    slong L = H->L;

    if (b - a <= L)
    {
        hyp_leaf(T, Q, R, FQ, FR, a, b, lev, need_r, 1, H);
    }
    else
    {
        mp_real_struct * T2 = H->tmp + 5 * lev;
        hyp_flist * f2 = H->gcd ? H->ftmp + 3 * lev : NULL;
        slong nl = (b - a + L - 1) / L, m = a + ((nl + 1) / 2) * L;

        int par = _hyp_par(nl, b - a, H);

        _hyp_halves(T, Q, R, FQ, FR, T2, f2, a, m, b, lev, need_r, 0,
            par == 1, H);
        if (H->gcd)
        {
            hyp_pending(H->ts->pl, lev, (m - 1) / L, H->ts);
            if (!hyp_remove(R, T2 + 1, FR, f2, H->ts->pl, f2 + 2,
                    H->J - lev, b - a, H->si, par == 2))
                hyp_stop(FQ, FR, f2, lev, H->ts);
        }
        (par == 2 ? merge_par : merge)(T, Q, R, T2, T2 + 1, T2 + 2,
            T2 + 3, T2 + 4, NULL, NULL, need_r, H->wp);
        if (H->gcd)
            hyp_lists_merge(FQ, FR, f2, f2 + 1, H->ts->pl, Q, R, need_r,
                &H->ts->scr);
        _hyp_release(T2, f2);
    }
}

/* ---- the tail bound ---- */

/* |a| / |b| (b != 0) as a double rounded up (up = 1) or down, for
   integers of any size: an upper bound is at least 2^-1000 |a|/|b|
   scaled into range (and HUGE_VAL beyond 2^1000), a lower bound 0 below
   2^-1000 */
static double
_ratio_bound(const fmpz_t a, const fmpz_t b, int up)
{
    slong ea, eb, e;
    double ma, mb, r;

    if (fmpz_is_zero(a))
        return 0.0;

    /* small integers convert exactly up to 2^-53 */
    if (!COEFF_IS_MPZ(*a) && !COEFF_IS_MPZ(*b))
    {
        r = (double) FLINT_UABS(*a) / (double) FLINT_UABS(*b);
        return up ? r * (1.0 + 0x1p-50) : r * (1.0 - 0x1p-50);
    }

    ma = fabs(fmpz_get_d_2exp(&ea, a));
    mb = fabs(fmpz_get_d_2exp(&eb, b));
    r = up ? ma / mb * (1.0 + 0x1p-48) : ma / mb * (1.0 - 0x1p-48);
    e = ea - eb;

    if (e > 1000)
        return up ? HUGE_VAL : d_mul_2exp_inrange(r, 1000);
    if (e < -1000)
        return up ? d_mul_2exp_inrange(r, -1000) : 0.0;
    return d_mul_2exp_inrange(r, (int) e);
}

/* a / b with the sign, as an estimate within the double range */
static double
_ratio_d(const fmpz_t a, const fmpz_t b)
{
    double r = _ratio_bound(a, b, 1);
    if (r == HUGE_VAL)
        r = 0x1p1000;
    return (fmpz_sgn(a) * fmpz_sgn(b) < 0) ? -r : r;
}

/* The tail after N terms is bounded as follows.  With every polynomial
   X of degree d written X(k) = lc k^d (1 + eps), |eps| <= e_X(k) =
   sum_{i<d} |c_i / lc| k^(i-d), decreasing in k, the ratios for k >= K
   satisfy |R(k)/Q(k)| <= rho(K), |P(k)/Q(k)| <= pi(K) (k/K)^delta,
   so that, with lp(N) an upper bound for log2 prod_{j<=N} |R(j)/Q(j)|,

       sum_{k>N} |term_k| <= 2^lp(N) pi(K) / (1 - rho(K) e^(delta/K)),

   K = N + 1, delta = max(0, deg P - deg Q).  lp(N) is accumulated term
   by term (the scan), from double Horner evaluations with a rigorous
   error bound (exact fmpz evaluations when that is inconclusive or a
   coefficient leaves the double range), up to the point K0 beyond
   which it is bounded in closed form (see _closed_log2); before the
   leading-term bounds apply, a bound from Taylor shifts serves (see
   tail_extra_shift).  The number of terms is found by the scan below
   K0 and by Newton steps on the closed form beyond, and the bound is
   sharp (to a fraction of a term) in both. */
typedef struct
{
    slong dP, dQ, dR;
    double rRQ, rPQ;            /* |lc(R)/lc(Q)|, |lc(P)/lc(Q)| (upper
                                   bounds) */
    double * eP, * eQ, * eR;    /* |c_i / lc| (upper bounds), i < d */
    const fmpz * Pf, * Qf, * Rf;
    slong Plen, Qlen, Rlen;
    double * Qd, * Rd, * Qa, * Ra;  /* coefficients, |coefficients|,
                                       scaled by a common 2^-sh */
    double * Pd;
    double geps;                /* Horner error, relative to the
                                   absolute evaluation */
    int dbl_ok;                 /* the scaled doubles are usable */
    /* the scan: prod_{j<=N} |R(j)/Q(j)| <= m 2^e */
    slong N;
    double m;
    slong e;
    int terminated;             /* R(j) = 0 for some j <= N */
    /* the closed form beyond K0 - 1 */
    double s1;                  /* upper bound for the 1/j coefficient
                                   of log |R(j)/Q(j)| (see below) */
    slong K0;                   /* WORD_MAX: none */
    double lpK;                 /* log2 prod_{j<K0} (upper bound) */
    double S2;                  /* second-order slack (nats) */
}
tail_struct;

static void
_rel_coeffs(double * c, slong * d, const fmpz * f, slong len)
{
    slong i;

    *d = len - 1;
    for (i = 0; i + 1 < len; i++)
        c[i] = _ratio_bound(f + i, f + len - 1, 1);
}

/* e_X(K) = sum_{i<d} c_i K^(i-d) (upper bound) */
static double
_rel_err(const double * c, slong d, double K)
{
    double v = 0.0, iK = 1.0 / K;
    slong i;

    for (i = 0; i < d; i++)
        v = (v + c[i]) * iK;
    return v * (1.0 + 1e-12);
}

/* The closed form.  For j >= K with e_Q(K), e_R(K) <= 1/2, writing
   R(j)/Q(j) = (r_d/q_d) j^(dR-dQ) (1 + u_R) / (1 + u_Q) with u_X(j) =
   sum_{i<d} (c_i/lc) j^(i-d),

       log |R(j)/Q(j)| <= log(lR/lQ) + (dR - dQ) log j + u_R - u_Q + u_Q^2

   (log(1+u) <= u, -log(1+u) <= -u + u^2 for |u| <= 1/2), and
   u_R - u_Q <= s1 / j + E2_R(j) + E2_Q(j), E2_X(j) = sum_{i<=d-2}
   |c_i/lc| j^(i-d) <= E2_X(K) (K/j)^2, e_Q(j)^2 <= e_Q(K)^2 (K/j)^2;
   the second-order terms sum to at most S2 = (E2_R(K) + E2_Q(K) +
   e_Q(K)^2) (K + 1) over j >= K.  K0 is the least such K (up to a
   factor two) with S2 below half a bit, and the product up to N >= K0
   - 1 is bounded by the exact scan to K0 - 1 plus the closed form. */

static double
_E2(const double * c, slong d, double K)
{
    return (d >= 2) ? _rel_err(c, d - 1, K) / K : 0.0;
}

static double
_S2(const tail_struct * t, double K)
{
    double eQ = _rel_err(t->eQ, t->dQ, K);
    return (_E2(t->eR, t->dR, K) + _E2(t->eQ, t->dQ, K) + eQ * eQ)
        * (K + 1.0) * (1.0 + 1e-12);
}

static int
_K0_ok(const tail_struct * t, double K)
{
    return _rel_err(t->eQ, t->dQ, K) <= 0.5
        && _rel_err(t->eR, t->dR, K) <= 0.5
        && _S2(t, K) <= 0.35;
}

static void
tail_choose_K0(tail_struct * t)
{
    slong K = 1, lo, hi, mid;

    while (!_K0_ok(t, (double) K))
    {
        if (K > (WORD(1) << (FLINT_BITS == 64 ? 40 : 29)))
        {
            t->K0 = WORD_MAX;
            return;
        }
        K *= 2;
    }

    /* the least K in (K/2, K] (the conditions are monotone) */
    lo = K / 2;
    hi = K;
    while (hi - lo > 1)
    {
        mid = lo + (hi - lo) / 2;
        if (_K0_ok(t, (double) mid))
            hi = mid;
        else
            lo = mid;
    }
    t->K0 = hi;
    t->S2 = _S2(t, (double) hi);
}

/* upper bound for log2 prod_{j=K0}^{N} |R(j)/Q(j)|, N >= K0 - 1 */
static double
_closed_log2(const tail_struct * t, slong N)
{
    double K = (double) t->K0, cnt = (double) (N - t->K0 + 1), v, g, ln2;

    if (N < t->K0)
        return 0.0;

    ln2 = 0.69314718055994530942;
    v = cnt * log(t->rRQ);
    v += fabs(v) * 1e-14;

    if (t->dR != t->dQ)
    {
        /* (dR - dQ) < 0 times a lower bound for sum_{j=K}^{N} log j */
        g = lgamma((double) N + 1.0) - lgamma(K);
        g = g * (1.0 - 1e-10) - 1e-9;
        v += (double) (t->dR - t->dQ) * FLINT_MAX(g, 0.0);
    }

    /* s1 sum_{j=K}^{N} 1/j */
    if (t->s1 >= 0.0)
        v += t->s1 * (1.0 / K + log((double) N / K) * (1.0 + 1e-14));
    else
        v += t->s1 * log(((double) N + 1.0) / K) * (1.0 - 1e-14);

    v += t->S2;
    return (v + fabs(v) * 1e-12) / ln2 + 1e-9;
}

static void
tail_init(tail_struct * t, const fmpz * P, slong Plen, const fmpz * Q,
    slong Qlen, const fmpz * R, slong Rlen)
{
    slong i;

    /* one allocation for all the double arrays */
    t->eP = flint_malloc((2 * Plen + 3 * Qlen + 3 * Rlen) * sizeof(double));
    t->eQ = t->eP + Plen;
    t->eR = t->eQ + Qlen;
    t->Qd = t->eR + Rlen;
    t->Rd = t->Qd + Qlen;
    t->Qa = t->Rd + Rlen;
    t->Ra = t->Qa + Qlen;
    t->Pd = t->Ra + Rlen;
    /* doubles scaled by a common 2^-sh (ratios unchanged), with the
       exact evaluations in tail_step as the fallback when a coefficient
       leaves the double range */
    {
        slong sh = 0, e;
        double m;

        for (i = 0; i < Plen; i++)
            sh = FLINT_MAX(sh, (slong) fmpz_bits(P + i));
        for (i = 0; i < Qlen; i++)
            sh = FLINT_MAX(sh, (slong) fmpz_bits(Q + i));
        for (i = 0; i < Rlen; i++)
            sh = FLINT_MAX(sh, (slong) fmpz_bits(R + i));
        sh = FLINT_MAX(sh - 900, 0);
        t->dbl_ok = 1;

#define SCALED(dst, x)                                              \
        do {                                                        \
            if (sh == 0 && !COEFF_IS_MPZ(*(x)))                     \
            {                                                       \
                (dst) = (double) *(x);                              \
                break;                                              \
            }                                                       \
            m = fmpz_get_d_2exp(&e, (x));                           \
            if (!fmpz_is_zero(x) && e - sh < -900)                  \
                t->dbl_ok = 0;                                      \
            (dst) = (e - sh < -1000) ? 0.0 : d_mul_2exp(m, (int) (e - sh)); \
        } while (0)

        for (i = 0; i < Plen; i++)
            SCALED(t->Pd[i], P + i);
        for (i = 0; i < Qlen; i++)
        {
            SCALED(t->Qd[i], Q + i);
            t->Qa[i] = fabs(t->Qd[i]);
        }
        for (i = 0; i < Rlen; i++)
        {
            SCALED(t->Rd[i], R + i);
            t->Ra[i] = fabs(t->Rd[i]);
        }
#undef SCALED
    }
    _rel_coeffs(t->eP, &t->dP, P, Plen);
    _rel_coeffs(t->eQ, &t->dQ, Q, Qlen);
    t->rPQ = _ratio_bound(P + Plen - 1, Q + Qlen - 1, 1);
    t->Pf = P;
    t->Plen = Plen;
    t->Qf = Q;
    t->Rf = R;
    t->Qlen = Qlen;
    t->Rlen = Rlen;
    t->N = 0;
    t->m = 1.0;
    t->e = 0;
    t->terminated = (Rlen == 0);

    t->K0 = WORD_MAX;

    if (Rlen == 0)
        return;

    _rel_coeffs(t->eR, &t->dR, R, Rlen);
    t->rRQ = _ratio_bound(R + Rlen - 1, Q + Qlen - 1, 1);

    if (t->dR > t->dQ || (t->dR == t->dQ && !(t->rRQ < 1.0)))
        flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: the series "
            "does not converge geometrically\n");

    /* coefficient conversion, Horner and the absolute evaluation */
    t->geps = (2.0 * FLINT_MAX(Qlen, Rlen) + 6.0) * 0x1p-52;

    /* s1 = r_{d-1} / r_d - q_{d-1} / q_d (the 1/j terms of the
       expansions), rounded up */
    {
        double a = (t->dR >= 1) ? _ratio_d(R + Rlen - 2, R + Rlen - 1)
            : 0.0;
        double b = (t->dQ >= 1) ? _ratio_d(Q + Qlen - 2, Q + Qlen - 1)
            : 0.0;
        t->s1 = (a - b) + (fabs(a) + fabs(b)) * 1e-14 + 1e-300;
    }
    tail_choose_K0(t);
}

static void
tail_clear(tail_struct * t)
{
    flint_free(t->eP);
}

/* (value, bound on |value - true|) of f(x) by Horner's rule */
static double
_horner_d(double * err, const double * f, const double * fa, slong len,
    double x, double geps)
{
    double v = 0.0, a = 0.0;
    slong i;

    for (i = len - 1; i >= 0; i--)
    {
        v = v * x + f[i];
        a = a * x + fa[i];
    }
    *err = a * geps;
    return v;
}

/* the scan: include the ratio at j = N + 1 */
static void
tail_step(tail_struct * t)
{
    slong j = t->N + 1;
    double r, q, er, eq, ratio;
    int ex;

    if (t->terminated)
    {
        t->N = j;
        return;
    }

    r = _horner_d(&er, t->Rd, t->Ra, t->Rlen, (double) j, t->geps);
    q = _horner_d(&eq, t->Qd, t->Qa, t->Qlen, (double) j, t->geps);
    r = fabs(r);
    q = fabs(q);

    if (t->dbl_ok && r > er && q > 2.0 * eq && r + er < 1e300
        && q < 1e300 && r > 1e-290)
    {
        /* both in the normal range: the quotient is too, and its
           exponent goes to e */
        ratio = frexp((r + er) / (q - eq) * (1.0 + 0x1p-50), &ex);
        t->e += ex;
    }
    else
    {
        fmpz_t rr, qq, kk;
        slong erx, eqx;
        double mr, mq;

        fmpz_init(rr);
        fmpz_init(qq);
        fmpz_init_set_ui(kk, j);
        _fmpz_poly_evaluate_fmpz(rr, t->Rf, t->Rlen, kk);
        _fmpz_poly_evaluate_fmpz(qq, t->Qf, t->Qlen, kk);
        if (fmpz_is_zero(qq))
            flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: Q(k) = 0\n");
        if (fmpz_is_zero(rr))
        {
            t->terminated = 1;
            ratio = 0.0;
        }
        else
        {
            mr = fabs(fmpz_get_d_2exp(&erx, rr)) * (1.0 + 0x1p-45);
            mq = fabs(fmpz_get_d_2exp(&eqx, qq)) * (1.0 - 0x1p-45);
            t->e += erx - eqx;
            ratio = mr / mq;
        }
        fmpz_clear(rr);
        fmpz_clear(qq);
        fmpz_clear(kk);
    }

    /* ratio in [1/4, 4] (or 0) */
    t->m *= ratio * (1.0 + 0x1p-50);
    if (!(t->m >= 0x1p-600 && t->m <= 0x1p600))
    {
        t->m = frexp(t->m, &ex);
        t->e += ex;
    }
    t->N = j;
}

/* upper bound for log2 sum_{k > N} |term_k| at the scan position N;
   HUGE_VAL if the leading-term bounds do not apply yet */
/* K^d for integer d (|d| small) */
static double
_pow_si(double K, slong d)
{
    double r = 1.0;
    slong i;
    for (i = 0; i < FLINT_ABS(d); i++)
        r *= K;
    return (d >= 0) ? r : 1.0 / r;
}

/* log2 (pi(K) / (1 - rho(K) e^(delta/K))), HUGE_VAL when the
   leading-term bounds do not apply at K */
static double
tail_extra(const tail_struct * t, double K)
{
    double rho, pi, eP, eQ, eR, b;
    slong delta = t->dP - t->dQ;

    eQ = _rel_err(t->eQ, t->dQ, K);
    if (!(eQ < 0.5))
        return HUGE_VAL;
    eR = _rel_err(t->eR, t->dR, K);
    eP = _rel_err(t->eP, t->dP, K);

    rho = (1.0 + eR) / (1.0 - eQ) * (t->rRQ)
        * _pow_si(K, t->dR - t->dQ);
    pi = (1.0 + eP) / (1.0 - eQ) * t->rPQ * _pow_si(K, delta);
    if (delta > 0)
        rho *= exp((double) delta / K);
    rho *= 1.000001;
    if (!(rho < 1.0))
        return HUGE_VAL;

    b = log2(pi / (1.0 - rho));
    return b + fabs(b) * 1e-12 + 1e-6;
}

/* The shifted bound, valid from small K when the leading-term bounds
   are not (Q with a large constant term, say).  With the exact Taylor
   shifts X(K + t) = sum_i x_i t^i, if the b_i = q_i all have the same
   sign then |Q(K + t)| = sum |b_i| t^i for t >= 0, whence

       |R(K + t)| <= rho sum |b_i| t^i,  rho = max_i |r_i| / |b_i|,
       |P(K + t)| <= pi (1 + t)^delta |Q(K + t)|,
                     pi = max_j (sum_{min(i, dQ) = j} |p_i|) / |b_j|,

   and the terms beyond N = K - 1 sum to at most 2^lp(N) pi delta! /
   (1 - rho)^(delta + 1) (as (1 + m)^delta <= delta! C(m + delta,
   delta)).  Returns the log2 of the factor after 2^lp(N), or
   HUGE_VAL. */
static double
tail_extra_shift(const tail_struct * t, slong K)
{
    slong Plen = t->Plen, Qlen = t->Qlen, Rlen = t->Rlen, i, j;
    slong dQ = Qlen - 1, delta = FLINT_MAX(Plen - Qlen, 0);
    fmpz * v;
    fmpz_t c;
    double rho = 0.0, pi = 0.0, r = HUGE_VAL;
    int sgn;

    if (Rlen > Qlen)
        return HUGE_VAL;

    v = _fmpz_vec_init(Plen + Qlen + Rlen);
    fmpz_init_set_ui(c, K);
    _fmpz_vec_set(v, t->Pf, Plen);
    _fmpz_vec_set(v + Plen, t->Qf, Qlen);
    _fmpz_vec_set(v + Plen + Qlen, t->Rf, Rlen);
    _fmpz_poly_taylor_shift(v, c, Plen);
    _fmpz_poly_taylor_shift(v + Plen, c, Qlen);
    _fmpz_poly_taylor_shift(v + Plen + Qlen, c, Rlen);

    sgn = fmpz_sgn(v + Plen + dQ);
    for (i = 0; i < Qlen; i++)
        if (fmpz_sgn(v + Plen + i) != sgn)
            goto cleanup;

    for (i = 0; i < Rlen; i++)
        rho = FLINT_MAX(rho, _ratio_bound(v + Plen + Qlen + i, v + Plen + i, 1));
    if (!(rho < 1.0))
        goto cleanup;

    for (j = 0; j <= dQ; j++)
    {
        double acc = 0.0;
        for (i = j; i < Plen && (i == j || j == dQ); i++)
            acc += _ratio_bound(v + i, v + Plen + j, 1);
        pi = FLINT_MAX(pi, acc * (1.0 + 1e-15));
    }
    if (pi == HUGE_VAL)
        goto cleanup;

    /* log2 (pi delta! / (1 - rho)^(delta + 1)) */
    r = log2(pi) - (double) (delta + 1) * log2(1.0 - rho);
    for (i = 2; i <= delta; i++)
        r += log2((double) i);
    r += fabs(r) * 1e-12 + 1e-6;

cleanup:
    _fmpz_vec_clear(v, Plen + Qlen + Rlen);
    fmpz_clear(c);
    return r;
}

/* upper bound for log2 sum_{k > N} |term_k| at the scan position N,
   with the shifted bound as a fallback when allowed */
static double
tail_log2(const tail_struct * t, int allow_shift)
{
    double x;

    if (t->terminated)
        return -HUGE_VAL;
    x = tail_extra(t, (double) (t->N + 1));
    if (x == HUGE_VAL && allow_shift)
        x = tail_extra_shift(t, t->N + 1);
    if (x == HUGE_VAL)
        return HUGE_VAL;
    return x + log2(t->m) + (double) t->e;
}

/* the same for N >= K0 - 1, by the closed form */
static double
tail_log2_closed(const tail_struct * t, slong N)
{
    double x = tail_extra(t, (double) (N + 1));
    if (x == HUGE_VAL)
        return HUGE_VAL;
    return x + t->lpK + _closed_log2(t, N);
}

/* the tail bound after N terms, N >= the scan position (advancing the
   scan when N < K0 - 1) */
static double
tail_bound(tail_struct * t, slong N)
{
    if (N < t->K0 - 1 || t->terminated)
    {
        while (t->N < N && !t->terminated)
            tail_step(t);
        return tail_log2(t, 1);
    }
    return tail_log2_closed(t, N);
}

/* 2^a for the mantissa thresholds of tail_find, whose mantissas stay
   within [2^-600, 2^600]: beyond +-1000 only the comparison matters,
   and no subnormal arises */
static double
_exp2_thr(double a)
{
    if (a < -1000.0)
        return 0.0;
    if (a > 1000.0)
        return HUGE_VAL;
    return exp2(a);
}

/* the least N (up to a small excess) with tail bound <= target, by the
   exact scan below K0 - 1 and Newton steps on the closed form beyond;
   returns N and sets *bound */
static slong
tail_find(double * bound, tail_struct * t, double target)
{
    double b, extra = HUGE_VAL, thr = 0.0, f, slope, cap;
    slong next = 0, next_shift = 16, e = WORD_MIN, N, Nprev, i;

    /* a generous cap on the term count: 4 times what the asymptotic
       rate needs, plus 2^20 (beyond it, the bounds cannot be
       established, e.g. for a Q with a huge positive root) */
    cap = 4.0 * (-target + 64.0);
    if (t->Rlen != 0 && t->dR == t->dQ)
        cap /= -log2(t->rRQ);
    cap = FLINT_MIN(cap + 1048576.0, 1e15);

    /* the scan */
    while (t->N < t->K0 - 1 || t->terminated)
    {
        if (t->e != e && extra < HUGE_VAL)
        {
            e = t->e;
            thr = _exp2_thr(target + 0.5 - extra - (double) e);
        }

        if (t->terminated || t->N >= next || t->m <= thr)
        {
            b = tail_log2(t, 0);
            /* the (costlier) shifted bound at geometric checkpoints and
               when predicted to succeed */
            if (b == HUGE_VAL && (t->N >= next_shift || t->m <= thr))
            {
                b = tail_log2(t, 1);
                next_shift = t->N + t->N / 4 + 1;
            }
            if (b <= target)
            {
                *bound = b;
                return t->N;
            }
            if (b < HUGE_VAL)
            {
                extra = b - (log2(t->m) + (double) t->e);
                e = t->e;
                thr = _exp2_thr(target + 0.5 - extra - (double) e);
                next = 2 * t->N + 1;
            }
            else
            {
                /* cheap: tail_extra fails early */
                next = t->N + 1;
            }
        }
        if ((double) t->N > cap)
            flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: unable to "
                "bound the tail\n");
        tail_step(t);
    }

    /* the closed form from N = K0 - 1 */
    t->lpK = log2(t->m) + (double) t->e;
    t->lpK += fabs(t->lpK) * 1e-12 + 1e-9;

    N = t->K0 - 1;
    Nprev = N;
    for (i = 0; ; i++)
    {
        b = tail_log2_closed(t, N);
        if (b <= target)
            break;
        if ((double) N > cap)
            flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: unable to "
                "bound the tail\n");
        Nprev = N;
        /* Newton step on the slope of the product at N */
        slope = log2(t->rRQ)
            + (double) (t->dR - t->dQ) * log2((double) N + 1.0);
        f = (b == HUGE_VAL) ? HUGE_VAL : b - target;
        if (slope < -1e-3 && f < 1e15)
            N += FLINT_MAX(1, (slong) ceil(f / -slope));
        else
            N += N + 1;
    }

    /* back off an overshoot in (Nprev, N]: probe below N by the
       predicted margin (usually just N - 1), bisecting when that
       prediction fails */
    {
        slong lo = Nprev, hi = N, c;
        double bh = b, bc;

        slope = log2(t->rRQ)
            + (double) (t->dR - t->dQ) * log2((double) N + 1.0);

        while (hi - lo > 1)
        {
            if (slope < -1e-3 && i >= 0)
            {
                c = hi - FLINT_MAX(1, (slong) floor((target - bh) / -slope));
                i = -1;     /* predict once */
            }
            else
                c = lo + (hi - lo) / 2;
            if (c <= lo)
                c = lo + (hi - lo) / 2;

            bc = tail_log2_closed(t, c);
            if (bc <= target)
            {
                hi = c;
                bh = bc;
            }
            else
                lo = c;
        }
        N = hi;
        b = bh;
    }

    *bound = b;
    return N;
}

/* ---- the driver ---- */

static void
_fmpz_set_signed_mpn(fmpz_t f, const mp_real_hypgeom_int_struct * x)
{
    if (x->n <= 0)
    {
        fmpz_zero(f);
        return;
    }
    fmpz_set_ui_array(f, x->d, x->n);
    if (x->neg)
        fmpz_neg(f, f);
}

/* x += [-2^e, 2^e] */
static void
_mp_real_add_error_2exp(mp_real_t x, slong e)
{
    slong q = e >> (FLINT_BITS == 64 ? 6 : 5);
    _mp_real_add_error_ulps_at(x, (double) (UWORD(1) << (e - q * FLINT_BITS)), q);
}

/* the content of a polynomial (positive gcd; 1 for the zero
   polynomial) and the quotient */
static void
_content(fmpz_t c, fmpz * q, const fmpz * f, slong len)
{
    _fmpz_vec_content(c, f, len);
    if (fmpz_is_zero(c))
        fmpz_one(c);
    _fmpz_vec_scalar_divexact_fmpz(q, f, len, c);
}

/* A = P / R exactly (Plen >= Rlen >= 1, lc(R) != 0); returns 0 if R
   does not divide P */
static int
_poly_divexact(fmpz * A, const fmpz * P, slong Plen, const fmpz * R,
    slong Rlen)
{
    slong i, j, Alen = Plen - Rlen + 1;
    fmpz * r;
    int ok = 1;

    r = _fmpz_vec_init(Plen);
    _fmpz_vec_set(r, P, Plen);

    for (i = Alen - 1; i >= 0 && ok; i--)
    {
        if (!fmpz_divisible(r + i + Rlen - 1, R + Rlen - 1))
        {
            ok = 0;
            break;
        }
        fmpz_divexact(A + i, r + i + Rlen - 1, R + Rlen - 1);
        for (j = 0; j < Rlen - 1; j++)
            fmpz_submul(r + i + j, A + i, R + j);
    }

    for (i = 0; i < Rlen - 1 && ok; i++)
        if (!fmpz_is_zero(r + i))
            ok = 0;

    _fmpz_vec_clear(r, Plen);
    return ok;
}

/* res = c x, c an integer (res may not alias x) */
static void
_mp_real_mul_fmpz(mp_real_t res, const mp_real_t x, const fmpz_t c, slong wp)
{
    if (fmpz_is_zero(c))
    {
        mp_real_zero(res);
    }
    else if (fmpz_size(c) == 1)
    {
        ulong u = COEFF_IS_MPZ(*c) ? COEFF_TO_PTR(*c)->_mp_d[0]
            : FLINT_UABS(*c);
        if (u == 1)
            mp_real_set(res, x);
        else
            mp_real_mul_ui(res, x, u, wp);
        if (fmpz_sgn(c) < 0)
            mp_real_neg(res, res);
    }
    else
    {
        mp_real_t t;
        mp_real_init(t);
        mp_real_set_fmpz(t, c);
        mp_real_mul(res, x, t, wp);
        mp_real_clear(t);
    }
}

/* |f(k)| by Horner's rule in doubles (an estimate), and an upper
   estimate of log2 sum |f_i| k^i */
static double
_eval_d(const double * f, slong len, double k)
{
    double v = 0.0;
    slong i;
    for (i = len - 1; i >= 0; i--)
        v = v * k + f[i];
    return v;
}

void
mp_real_hypgeom_series(mp_real_t res, const mp_real_hypgeom_series_struct * s,
    slong n)
{
    _mp_real_hypgeom_series(res, s, n, -1);
}

void
_mp_real_hypgeom_series(mp_real_t res, const mp_real_hypgeom_series_struct * s,
    slong n, int gcd)
{
    hyp_struct H;
    hyp_sieve_info si;
    hyp_thread ts;
    hyp_flist FQt, FRt;
    tail_struct tl;
    fmpz * P, * Q, * R, * Qc, * Rc, * A = NULL;
    double * Pd, * Qd, * Rd;
    fmpz_t coefP, coefQ, coefD;
    slong Plen = s->Plen, Qlen = s->Qlen, Rlen = s->Rlen, Alen = 0;
    slong N, K, J, L, i, wp, depth, Lmax;
    double target, Sest, lS, tau, pt;
    mp_real_t T, Qt, Rt, num, den, t;

    if (s->power != 1 && s->power != -1)
        flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: power must be 1 or -1\n");

    wp = FLINT_MAX(n, 2);

    P = _fmpz_vec_init(Plen + Qlen + Rlen + 3 + Qlen + Rlen + 2);
    Q = P + Plen + 1;
    R = Q + Qlen + 1;
    Qc = R + Rlen + 1;
    Rc = Qc + Qlen + 1;
    fmpz_init(coefP);
    fmpz_init(coefQ);
    fmpz_init(coefD);
    fmpz_init(H.cQ);
    fmpz_init(H.cR);

    for (i = 0; i < Plen; i++)
        _fmpz_set_signed_mpn(P + i, s->P + i);
    for (i = 0; i < Qlen; i++)
        _fmpz_set_signed_mpn(Q + i, s->Q + i);
    for (i = 0; i < Rlen; i++)
        _fmpz_set_signed_mpn(R + i, s->R + i);
    _fmpz_set_signed_mpn(coefP, &s->coefP);
    _fmpz_set_signed_mpn(coefQ, &s->coefQ);
    _fmpz_set_signed_mpn(coefD, &s->coefD);

    while (Qlen > 0 && fmpz_is_zero(Q + Qlen - 1))
        Qlen--;
    if (Qlen == 0)
        flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: Q = 0\n");
    if (fmpz_is_zero(coefD))
        flint_throw(FLINT_ERROR, "mp_real_hypgeom_series: coefD = 0\n");
    while (Plen > 0 && fmpz_is_zero(P + Plen - 1))
        Plen--;
    while (Rlen > 0 && fmpz_is_zero(R + Rlen - 1))
        Rlen--;

    mp_real_init(t);

    /* the trivial sum */
    if (Plen == 0 || fmpz_is_zero(coefP))
    {
        mp_real_set_fmpz(res, coefQ);
        mp_real_set_fmpz(t, coefD);
        if (s->power == 1)
            mp_real_div(res, res, t, wp);
        else
            mp_real_div(res, t, res, wp);
        mp_real_clear(t);
        goto cleanup_early;
    }

    tail_init(&tl, P, Plen, Q, Qlen, R, Rlen);
    Pd = tl.Pd;
    Qd = tl.Qd;
    Rd = tl.Rd;

    /* magnitude estimate of S coefD / coefP from the first terms (only
       used to choose the target) */
    {
        double sum = 0.0, prod = 1.0, q;
        for (i = 1; i <= 40; i++)
        {
            q = _eval_d(Qd, Qlen, i);
            if (q == 0.0)
                break;
            sum += _eval_d(Pd, Plen, i) / q * prod;
            prod *= _eval_d(Rd, Rlen, i) / q;
            if (!(fabs(prod) > 1e-30 * fabs(sum)) || !isfinite(prod))
                break;
        }
        Sest = fabs(_ratio_d(coefQ, coefP) + sum);
        if (!(Sest > 0.0) || !isfinite(Sest))
            Sest = 1.0;
        if (fabs(_ratio_d(coefQ, coefP)) >= 0x1p1000)
            Sest = 0x1p1000;
        lS = log2(Sest);
        if (fabs(_ratio_d(coefQ, coefP)) >= 0x1p1000)
            lS = (double) fmpz_bits(coefQ) - (double) fmpz_bits(coefP);
    }

    /* target: tail below 2^-(FLINT_BITS wp + 16) Sest */
    target = lS - (double) (FLINT_BITS * wp + 16);

    /* the number of terms */
    N = tail_find(&tau, &tl, target);
    N = FLINT_MAX(N, 1);

    /* leaves of L terms: about HYP_LEAF_LIMBS limbs of T, between
       HYP_LEAF_TERMS_MIN and HYP_LEAF_TERMS_MAX terms; pt is an upper
       estimate of the bits per term */
    pt = FLINT_MAX(_log2_bound_fmpz(Q, Qlen, N),
        _log2_bound_fmpz(R, Rlen, N));
    pt = FLINT_MAX(pt, 1.0) + 2.0;
    Lmax = (slong) ((FLINT_BITS * HYP_LEAF_LIMBS) / pt);
    Lmax = FLINT_MAX(Lmax, HYP_LEAF_TERMS_MIN);
    Lmax = FLINT_MIN(Lmax, HYP_LEAF_TERMS_MAX);

    /* content mode */
    _content(H.cQ, Qc, Q, Qlen);
    _content(H.cR, Rc, R, Rlen);
    {
        double cbits = (double) fmpz_bits(H.cQ) - 0.5, qbits;
        qbits = _log2_bound_fmpz(Q, Qlen, N) - cbits;
        H.content = (N > Lmax && wp >= 64 && cbits >= 8.0 * qbits);
    }

    /* in content mode, the same rule for R: separate its content only
       when that dominates (otherwise R stays short against Q and T) */
    H.have_cR = 0;
    if (H.content)
    {
        double rc = (double) fmpz_bits(H.cR) - 0.5, rbits;
        rbits = _log2_bound_fmpz(R, Rlen, N) - rc;
        H.have_cR = !fmpz_is_one(H.cR) && rc >= 8.0 * rbits;
    }

    /* content removal: forced, or chosen by the precision and the share
       of the bits that can cancel; it makes the levels cheaper relative
       to the leaves, so that the leaves are shorter: about
       HYP_GCD_LEAF_LIMBS limbs, at least HYP_GCD_LEAF_TERMS_MIN terms
       (zeta(3), whose terms have five-limb factors, was measured 4%
       faster at 12 terms than at 24; Catalan, log 2 and pi about the
       same) */
    if (gcd < 0)
        gcd = (wp >= _hyp_gcd_min_limbs(1.0) && N > Lmax
            && wp >= _hyp_gcd_min_limbs(_hyp_gcd_share(Q, Qlen, R, Rlen, N,
                !H.content, !H.have_cR)));
    if (gcd && !H.content)
    {
        Lmax = (slong) ((FLINT_BITS * HYP_GCD_LEAF_LIMBS) / pt);
        Lmax = FLINT_MAX(Lmax, HYP_GCD_LEAF_TERMS_MIN);
        Lmax = FLINT_MIN(Lmax, HYP_LEAF_TERMS_MAX);
    }

    if (H.content)
    {
        /* a balanced tree of 2^J leaves, N rounded up to K L */
        for (J = 0, K = 1; K * Lmax < N; J++, K *= 2)
            ;
        L = (N + K - 1) / K;
        N = K * L;
    }
    else
    {
        L = Lmax;
        K = (N + L - 1) / L;
        J = (K <= 1) ? 0 : FLINT_BIT_COUNT(K - 1);
    }

    H.gcd = 0;
    H.si = NULL;
    H.ts = NULL;
    H.ftmp = NULL;
    if (J >= 1 && gcd)
        H.gcd = sieve_info_init(&si, Q, Qlen, R, Rlen, !H.content,
            !H.have_cR, N, L);

    /* factored mode: P = A R */
    H.factored = 0;
    if (Rlen > 0 && Plen >= Rlen)
    {
        Alen = Plen - Rlen + 1;
        A = _fmpz_vec_init(Alen);
        H.factored = _poly_divexact(A, P, Plen, R, Rlen);
    }

    hpoly_init(H.P, P, Plen, N);
    hpoly_init(H.Q, Q, Qlen, N);
    hpoly_init(H.R, R, Rlen, N);
    hpoly_init(H.A, A, H.factored ? Alen : 0, N);
    hpoly_init(H.Qc, Qc, H.content ? Qlen : 0, N);
    hpoly_init(H.Rc, Rc, H.content ? Rlen : 0, N);

    /* buffers: numbers of at most L terms, with room for the product
       by a value of P or A */
    pt = FLINT_MAX(pt, FLINT_BITS * FLINT_MAX(H.Q->W, H.R->W));
    H.bufn = (slong) ((L + 1) * pt) / FLINT_BITS
        + H.P->W + H.A->W + H.Q->W + H.R->W + 8;
    H.vlimbs = H.P->W + 2 * H.Q->W + 2 * H.R->W + H.A->W + 8;

    /* the tail at the final N (rounded up in content mode) */
    tau = FLINT_MIN(tau, tail_bound(&tl, N));

    H.L = L;
    H.J = J;
    H.wp = wp;
    H.tbits = FLINT_MAX(FLINT_MAX(_log2_bound_fmpz(Q, Qlen, N),
        _log2_bound_fmpz(R, Rlen, N)), 1.0) + 2.0;
    if (H.content)
        H.tbits = FLINT_MAX(_log2_bound_fmpz(Q, Qlen, N)
            - ((double) fmpz_bits(H.cQ) - 1.0), 1.0) + 2.0;

    mp_real_init(T);
    mp_real_init(Qt);
    mp_real_init(Rt);
    mp_real_init(num);
    mp_real_init(den);

    /* one allocation: the leaf buffers, then the coefficients */
    {
        slong csz = H.P->len * H.P->W + H.Q->len * H.Q->W
            + H.R->len * H.R->W + H.A->len * H.A->W
            + H.Qc->len * H.Qc->W + H.Rc->len * H.Rc->W;
        nn_ptr c;

        H.buf = flint_malloc((NBUF * H.bufn + H.vlimbs + csz)
            * sizeof(ulong));
        c = H.buf + NBUF * H.bufn + H.vlimbs;
        hpoly_fill(H.P, P, c); c += H.P->len * H.P->W;
        hpoly_fill(H.Q, Q, c); c += H.Q->len * H.Q->W;
        hpoly_fill(H.R, R, c); c += H.R->len * H.R->W;
        hpoly_fill(H.A, A, c); c += H.A->len * H.A->W;
        hpoly_fill(H.Qc, Qc, c); c += H.Qc->len * H.Qc->W;
        hpoly_fill(H.Rc, Rc, c);
    }

    if (J == 0)
    {
        /* a single exact leaf */
        leaf(T, Qt, Rt, 1, N + 1, 0, 1, &H);
    }
    else
    {
        depth = J + 1;
        H.tmp = flint_malloc(5 * depth * sizeof(mp_real_struct));
        for (i = 0; i < 5 * depth; i++)
            mp_real_init(H.tmp + i);
        flist_init(&FQt);
        flist_init(&FRt);
        if (H.gcd)
        {
            H.si = &si;
            thread_init(&ts, depth);
            H.ts = &ts;
            H.ftmp = flint_malloc(3 * depth * sizeof(hyp_flist));
            for (i = 0; i < 3 * depth; i++)
                flist_init(H.ftmp + i);
        }

        if (H.content)
        {
            fmpz_t c;
            fmpz_init(c);
            H.cQpow = flint_malloc(2 * J * sizeof(mp_real_struct));
            H.cRpow = H.cQpow + J;
            for (i = 0; i < 2 * J; i++)
                mp_real_init(H.cQpow + i);
            fmpz_pow_ui(c, H.cQ, L);
            mp_real_set_fmpz(H.cQpow, c);
            for (i = 1; i < J; i++)
                mp_real_mul(H.cQpow + i, H.cQpow + i - 1, H.cQpow + i - 1, wp);
            if (H.have_cR)
            {
                fmpz_pow_ui(c, H.cR, L);
                mp_real_set_fmpz(H.cRpow, c);
                for (i = 1; i < J; i++)
                    mp_real_mul(H.cRpow + i, H.cRpow + i - 1, H.cRpow + i - 1, wp);
            }
            fmpz_clear(c);

            bsplit_pow2(T, Qt, Rt, &FQt, &FRt, 0, K, J, 0, &H);

            /* Q = Q' cQ^N, cQ^N = (cQ^(N/2))^2 */
            mp_real_mul(t, H.cQpow + (J - 1), H.cQpow + (J - 1), wp);
            mp_real_mul(num, Qt, t, wp);
            mp_real_swap(Qt, num);
            for (i = 0; i < 2 * J; i++)
                mp_real_clear(H.cQpow + i);
            flint_free(H.cQpow);
        }
        else
        {
            bsplit_range(T, Qt, Rt, &FQt, &FRt, 1, N + 1, J, 0, &H);
        }

        for (i = 0; i < 5 * depth; i++)
            mp_real_clear(H.tmp + i);
        flint_free(H.tmp);
        flist_clear(&FQt);
        flist_clear(&FRt);
        if (H.gcd)
        {
            for (i = 0; i < 3 * depth; i++)
                flist_clear(H.ftmp + i);
            flint_free(H.ftmp);
            thread_clear(&ts);
        }
    }
    if (H.gcd)
        sieve_info_clear(&si);

    flint_free(H.buf);

    /* num = coefQ Q + coefP T (+/- |coefP| 2^tau |Q|), den = coefD Q */
    _mp_real_mul_fmpz(num, T, coefP, wp);
    if (tau > -HUGE_VAL)
    {
        /* |Q| < B^exp */
        double e = tau + (double) fmpz_bits(coefP) + 1.0;
        _mp_real_add_error_2exp(num, (slong) ceil(e) + FLINT_BITS * Qt->exp);
    }
    if (!fmpz_is_zero(coefQ))
    {
        _mp_real_mul_fmpz(T, Qt, coefQ, wp);
        mp_real_add(num, num, T, wp);
    }

    if (fmpz_is_one(coefD))
    {
        mp_real_swap(den, Qt);
    }
    else
    {
        _mp_real_mul_fmpz(den, Qt, coefD, wp);
    }

    if (s->power == 1)
        mp_real_div(res, num, den, wp);
    else
        mp_real_div(res, den, num, wp);

    mp_real_clear(T);
    mp_real_clear(Qt);
    mp_real_clear(Rt);
    mp_real_clear(num);
    mp_real_clear(den);
    mp_real_clear(t);

    tail_clear(&tl);
    if (A != NULL)
        _fmpz_vec_clear(A, Alen);

cleanup_early:
    _fmpz_vec_clear(P, s->Plen + s->Qlen + s->Rlen + 3 + s->Qlen + s->Rlen + 2);
    fmpz_clear(coefP);
    fmpz_clear(coefQ);
    fmpz_clear(coefD);
    fmpz_clear(H.cQ);
    fmpz_clear(H.cR);
}

/* ---- the int64 convenience wrapper ---- */

/* x = v as a signed mpn integer in d (room for 64 / FLINT_BITS limbs) */
static void
_int_set_int64(mp_real_hypgeom_int_struct * x, nn_ptr d, int64_t v)
{
    uint64_t u = (v < 0) ? -(uint64_t) v : (uint64_t) v;

    x->neg = (v < 0);
    x->d = d;
#if FLINT_BITS == 64
    d[0] = u;
    x->n = (u != 0);
#else
    d[0] = (ulong) u;
    d[1] = (ulong) (u >> 32);
    x->n = (d[1] != 0) ? 2 : (d[0] != 0);
#endif
}

#define I64_LIMBS (64 / FLINT_BITS)

void
mp_real_hypgeom_series_int64(mp_real_t res, int power, int64_t coefP,
    int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen,
    const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n)
{
    mp_real_hypgeom_series_struct s;
    mp_real_hypgeom_int_struct * c;
    nn_ptr d;
    slong i, tot = Plen + Qlen + Rlen + 3;

    c = flint_malloc(tot * sizeof(mp_real_hypgeom_int_struct));
    d = flint_malloc(tot * I64_LIMBS * sizeof(ulong));

    for (i = 0; i < Plen; i++)
        _int_set_int64(c + i, d + i * I64_LIMBS, P[i]);
    for (i = 0; i < Qlen; i++)
        _int_set_int64(c + Plen + i, d + (Plen + i) * I64_LIMBS, Q[i]);
    for (i = 0; i < Rlen; i++)
        _int_set_int64(c + Plen + Qlen + i,
            d + (Plen + Qlen + i) * I64_LIMBS, R[i]);
    _int_set_int64(&s.coefP, d + (tot - 3) * I64_LIMBS, coefP);
    _int_set_int64(&s.coefQ, d + (tot - 2) * I64_LIMBS, coefQ);
    _int_set_int64(&s.coefD, d + (tot - 1) * I64_LIMBS, coefD);

    s.power = power;
    s.P = c;
    s.Plen = Plen;
    s.Q = c + Plen;
    s.Qlen = Qlen;
    s.R = c + Plen + Qlen;
    s.Rlen = Rlen;

    mp_real_hypgeom_series(res, &s, n);

    flint_free(c);
    flint_free(d);
}
