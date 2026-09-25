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
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fixed.h"

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
   final quotient is a single fball division, with coefP, coefQ, coefD
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
    fball_struct * cQpow;   /* cQ^(L 2^i), i < J */
    fball_struct * cRpow;
    fball_struct * tmp;     /* 5 per level */
    double tbits;           /* estimated bits per term of T, Q */
}
hyp_struct;

#define NBUF 11

/* the leaf when every value is a single word (and no content powers):
   the common case, with the per-term work reduced to the mpn passes */
static void
leaf_w1(fball_t T, fball_t Qo, fball_t Ro, slong a, slong b, int need_r,
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

    fball_set_mpn_2exp(T, Tb, tn, 0);
    if (tneg)
        fball_neg(T);
    if (content)
    {
        fball_set_mpn_2exp(Qo, Cs, cn, 0);
        if (cneg)
            fball_neg(Qo);
    }
    else
    {
        fball_set_mpn_2exp(Qo, Qs, qn, 0);
        if (qneg)
            fball_neg(Qo);
    }
    if (need_r)
    {
        fball_set_mpn_2exp(Ro, Rs, rn, 0);
        if (rneg)
            fball_neg(Ro);
    }
}

/* exact T, Q (or Q'), R (or R') over [a, b) as fballs; single-word
   polynomial values are computed for the whole block up front */
static void
leaf(fball_t T, fball_t Qo, fball_t Ro, slong a, slong b, int need_r,
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

    fball_set_mpn_2exp(T, Tb, tn, 0);
    if (tneg)
        fball_neg(T);

    if (content)
    {
        fball_set_mpn_2exp(Qo, Pa, pn, 0);
        if (pneg)
            fball_neg(Qo);
    }
    else
    {
        fball_set_mpn_2exp(Qo, Qa, qn, 0);
        if (qneg)
            fball_neg(Qo);
    }

    if (need_r)
    {
        fball_set_mpn_2exp(Ro, Ra, rn, 0);
        if (rneg)
            fball_neg(Ro);
    }
}

/* T = T1 Q2 + R1 T2 (with the content powers pq, pr when non-NULL),
   Q = Q1 Q2, R = R1 R2; (T, Q, R) hold the left values on entry;
   (T2, Q2, R2, U, V) is scratch */
static void
merge(fball_t T, fball_t Q, fball_t R, fball_t T2, fball_t Q2, fball_t R2,
    fball_t U, fball_t V, const fball_struct * pq, const fball_struct * pr,
    int need_r, slong wp)
{
    fball_mul(U, T, Q2, wp);
    if (pq != NULL)
    {
        fball_mul(V, U, pq, wp);
        fball_swap(U, V);
    }

    fball_mul(V, R, T2, wp);
    if (pr != NULL)
    {
        fball_mul(T2, V, pr, wp);
        fball_swap(T2, V);
    }

    fball_add(T, U, V, wp);

    fball_mul(U, Q, Q2, wp);
    fball_swap(Q, U);

    if (need_r)
    {
        fball_mul(U, R, R2, wp);
        fball_swap(R, U);
    }
}

static void bsplit_pow2(fball_t T, fball_t Q, fball_t R, slong a0, slong k,
    slong lev, int need_r, const hyp_struct * H);
static void bsplit_range(fball_t T, fball_t Q, fball_t R, slong a, slong b,
    slong lev, int need_r, const hyp_struct * H);

/* a subtree as a job for _fixed_parallel_pair; with own set it runs on
   private leaf buffers and level temporaries (the shared H only
   supplies the read-only polynomials and content powers), freed at the
   end */
typedef struct
{
    fball_struct * T, * Q, * R;
    slong a, b, lev;
    int need_r, pow2, own;
    const hyp_struct * H;
}
hyp_job;

static void
_hyp_job(void * arg)
{
    hyp_job * J = (hyp_job *) arg;
    hyp_struct H2;
    const hyp_struct * H = J->H;
    slong i, nt = 5 * (J->lev + 1);

    if (J->own)
    {
        H2 = *J->H;
        H2.buf = flint_malloc((NBUF * H2.bufn + H2.vlimbs) * sizeof(ulong));
        H2.tmp = flint_malloc(nt * sizeof(fball_struct));
        for (i = 0; i < nt; i++)
            fball_init(H2.tmp + i);
        H = &H2;
    }

    if (J->pow2)
        bsplit_pow2(J->T, J->Q, J->R, J->a, J->b, J->lev, J->need_r, H);
    else
        bsplit_range(J->T, J->Q, J->R, J->a, J->b, J->lev, J->need_r, H);

    if (J->own)
    {
        for (i = 0; i < nt; i++)
            fball_clear(H2.tmp + i);
        flint_free(H2.tmp);
        flint_free(H2.buf);
    }
}

/* subtrees of at least HYP_PAR_LEAVES leaves may run on two threads,
   by the rule of FIXED_PAR_CAP: 0 serial, 1 the halves on two
   threads, 2 above the cap (the merge on two threads) */
#define HYP_PAR_LEAVES 64

static int
_hyp_par(slong leaves, slong terms, const hyp_struct * H)
{
    if (leaves < HYP_PAR_LEAVES || flint_get_num_threads() < 2)
        return 0;
    return ((double) terms * H->tbits
        > FIXED_PAR_CAP * FLINT_BITS * (double) H->wp) ? 2 : 1;
}

static void
_hyp_halves(fball_struct * T, fball_struct * Q, fball_struct * R_,
    fball_struct * T2, slong a, slong m, slong b,
    slong lev, int need_r, int pow2, int fork, const hyp_struct * H)
{
    if (fork)
    {
        hyp_job L, R;
        L.T = T; L.Q = Q; L.R = R_; L.a = a; L.b = m;
        L.lev = lev - 1; L.need_r = 1; L.pow2 = pow2; L.own = 0; L.H = H;
        R.T = T2; R.Q = T2 + 1; R.R = T2 + 2; R.a = pow2 ? a + m : m;
        R.b = pow2 ? m : b; R.lev = lev - 1; R.need_r = need_r;
        R.pow2 = pow2; R.own = 1; R.H = H;
        _fixed_parallel_pair(_hyp_job, &L, _hyp_job, &R);
    }
    else if (pow2)
    {
        bsplit_pow2(T, Q, R_, a, m, lev - 1, 1, H);
        bsplit_pow2(T2, T2 + 1, T2 + 2, a + m, m, lev - 1, need_r, H);
    }
    else
    {
        bsplit_range(T, Q, R_, a, m, lev - 1, 1, H);
        bsplit_range(T2, T2 + 1, T2 + 2, m, b, lev - 1, need_r, H);
    }
}

/* merge on two threads: {U = T1 Q2 (pq), Q = Q1 Q2} and {V = R1 T2 (pr),
   R = R1 R2} */
typedef struct
{
    fball_struct * T, * Q, * R, * T2, * Q2, * R2, * U, * V;
    const fball_struct * pq, * pr;
    int need_r;
    slong wp;
}
hyp_merge_struct;

static void
_hyp_merge_x(void * arg)
{
    hyp_merge_struct * M = (hyp_merge_struct *) arg;
    fball_mul(M->U, M->T, M->Q2, M->wp);
    if (M->pq != NULL)
        fball_mul(M->U, M->U, M->pq, M->wp);
    fball_mul(M->Q, M->Q, M->Q2, M->wp);
}

static void
_hyp_merge_y(void * arg)
{
    hyp_merge_struct * M = (hyp_merge_struct *) arg;
    fball_mul(M->V, M->R, M->T2, M->wp);
    if (M->pr != NULL)
        fball_mul(M->V, M->V, M->pr, M->wp);
    if (M->need_r)
        fball_mul(M->R, M->R, M->R2, M->wp);
}

static void
merge_par(fball_t T, fball_t Q, fball_t R, fball_t T2, fball_t Q2,
    fball_t R2, fball_t U, fball_t V, const fball_struct * pq,
    const fball_struct * pr, int need_r, slong wp)
{
    hyp_merge_struct M;
    M.T = T; M.Q = Q; M.R = R; M.T2 = T2; M.Q2 = Q2; M.R2 = R2;
    M.U = U; M.V = V; M.pq = pq; M.pr = pr; M.need_r = need_r; M.wp = wp;
    _fixed_parallel_pair(_hyp_merge_x, &M, _hyp_merge_y, &M);
    fball_add(T, U, V, wp);
}

/* free this level's temporaries when large: they are otherwise held
   through all the merges above */
#define HYP_KEEP 2048

static void
_hyp_release(fball_struct * T2)
{
    slong i;
    if (T2[0].alloc > HYP_KEEP || T2[1].alloc > HYP_KEEP)
        for (i = 0; i < 5; i++)
        {
            fball_clear(T2 + i);
            fball_init(T2 + i);
        }
}

/* content mode: k leaves (a power of two) from leaf a0, at level
   lev = log2(k) */
static void
bsplit_pow2(fball_t T, fball_t Q, fball_t R, slong a0, slong k, slong lev,
    int need_r, const hyp_struct * H)
{
    if (k == 1)
    {
        leaf(T, Q, R, 1 + a0 * H->L, 1 + (a0 + 1) * H->L, need_r, 0, H);
    }
    else
    {
        fball_struct * T2 = H->tmp + 5 * lev;
        slong h = k / 2;

        int par = _hyp_par(k, k * H->L, H);

        _hyp_halves(T, Q, R, T2, a0, h, 0, lev, need_r, 1, par == 1, H);
        (par == 2 ? merge_par : merge)(T, Q, R, T2, T2 + 1, T2 + 2,
            T2 + 3, T2 + 4, H->cQpow + (lev - 1),
            H->have_cR ? H->cRpow + (lev - 1) : NULL, need_r, H->wp);
        _hyp_release(T2);
    }
}

/* carry mode: the terms [a, b), in leaves of L terms (the last one
   possibly shorter), at most 2^lev leaves */
static void
bsplit_range(fball_t T, fball_t Q, fball_t R, slong a, slong b, slong lev,
    int need_r, const hyp_struct * H)
{
    slong L = H->L;

    if (b - a <= L)
    {
        leaf(T, Q, R, a, b, need_r, 1, H);
    }
    else
    {
        fball_struct * T2 = H->tmp + 5 * lev;
        slong nl = (b - a + L - 1) / L, m = a + ((nl + 1) / 2) * L;

        int par = _hyp_par(nl, b - a, H);

        _hyp_halves(T, Q, R, T2, a, m, b, lev, need_r, 0, par == 1, H);
        (par == 2 ? merge_par : merge)(T, Q, R, T2, T2 + 1, T2 + 2,
            T2 + 3, T2 + 4, NULL, NULL, need_r, H->wp);
        _hyp_release(T2);
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
        return up ? HUGE_VAL : ldexp(r, 1000);
    if (e < -1000)
        return up ? ldexp(r, -1000) : 0.0;
    return ldexp(r, (int) e);
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
        if (K > (WORD(1) << 40))
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
            (dst) = (e - sh < -1000) ? 0.0 : ldexp(m, (int) (e - sh)); \
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
        flint_throw(FLINT_ERROR, "fball_hypgeom_series: the series "
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
            flint_throw(FLINT_ERROR, "fball_hypgeom_series: Q(k) = 0\n");
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
            thr = exp2(target + 0.5 - extra - (double) e);
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
                thr = exp2(target + 0.5 - extra - (double) e);
                next = 2 * t->N + 1;
            }
            else
            {
                /* cheap: tail_extra fails early */
                next = t->N + 1;
            }
        }
        if ((double) t->N > cap)
            flint_throw(FLINT_ERROR, "fball_hypgeom_series: unable to "
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
            flint_throw(FLINT_ERROR, "fball_hypgeom_series: unable to "
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
_fmpz_set_signed_mpn(fmpz_t f, const fixed_hypgeom_int_struct * x)
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

static void
_fball_set_fmpz(fball_t x, const fmpz_t f)
{
    slong n = fmpz_size(f);
    fmpz_t a;
    nn_ptr d;

    if (n == 0)
    {
        fball_zero(x);
        return;
    }
    fmpz_init(a);
    fmpz_abs(a, f);
    d = flint_malloc(n * sizeof(ulong));
    fmpz_get_ui_array(d, n, a);
    fball_set_mpn_2exp(x, d, n, 0);
    if (fmpz_sgn(f) < 0)
        fball_neg(x);
    flint_free(d);
    fmpz_clear(a);
}

/* x += [-2^e, 2^e] */
static void
_fball_add_error_2exp(fball_t x, slong e)
{
    slong q = e >> (FLINT_BITS == 64 ? 6 : 5);
    fball_add_error(x, ldexp(1.0, (int) (e - q * FLINT_BITS)), q);
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
_fball_mul_fmpz(fball_t res, const fball_t x, const fmpz_t c, slong wp)
{
    if (fmpz_is_zero(c))
    {
        fball_zero(res);
    }
    else if (fmpz_size(c) == 1)
    {
        ulong u = COEFF_IS_MPZ(*c) ? COEFF_TO_PTR(*c)->_mp_d[0]
            : FLINT_UABS(*c);
        if (u == 1)
            fball_set(res, x);
        else
            fball_mul_ui(res, x, u, wp);
        if (fmpz_sgn(c) < 0)
            fball_neg(res);
    }
    else
    {
        fball_t t;
        fball_init(t);
        _fball_set_fmpz(t, c);
        fball_mul(res, x, t, wp);
        fball_clear(t);
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
fball_hypgeom_series(fball_t res, const fixed_hypgeom_series_struct * s,
    slong n)
{
    hyp_struct H;
    tail_struct tl;
    fmpz * P, * Q, * R, * Qc, * Rc, * A = NULL;
    double * Pd, * Qd, * Rd;
    fmpz_t coefP, coefQ, coefD;
    slong Plen = s->Plen, Qlen = s->Qlen, Rlen = s->Rlen, Alen = 0;
    slong N, K, J, L, i, wp, depth, Lmax;
    double target, Sest, lS, tau, pt;
    fball_t T, Qt, Rt, num, den, t;

    if (s->power != 1 && s->power != -1)
        flint_throw(FLINT_ERROR, "fball_hypgeom_series: power must be 1 or -1\n");

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
        flint_throw(FLINT_ERROR, "fball_hypgeom_series: Q = 0\n");
    if (fmpz_is_zero(coefD))
        flint_throw(FLINT_ERROR, "fball_hypgeom_series: coefD = 0\n");
    while (Plen > 0 && fmpz_is_zero(P + Plen - 1))
        Plen--;
    while (Rlen > 0 && fmpz_is_zero(R + Rlen - 1))
        Rlen--;

    fball_init(t);

    /* the trivial sum */
    if (Plen == 0 || fmpz_is_zero(coefP))
    {
        _fball_set_fmpz(res, coefQ);
        _fball_set_fmpz(t, coefD);
        if (s->power == 1)
            fball_div(res, res, t, wp);
        else
            fball_div(res, t, res, wp);
        fball_clear(t);
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
            Sest = ldexp(1.0, 1000);
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

    H.have_cR = 0;

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

    fball_init(T);
    fball_init(Qt);
    fball_init(Rt);
    fball_init(num);
    fball_init(den);

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
        H.tmp = flint_malloc(5 * depth * sizeof(fball_struct));
        for (i = 0; i < 5 * depth; i++)
            fball_init(H.tmp + i);

        if (H.content)
        {
            fmpz_t c;
            fmpz_init(c);
            H.cQpow = flint_malloc(2 * J * sizeof(fball_struct));
            H.cRpow = H.cQpow + J;
            for (i = 0; i < 2 * J; i++)
                fball_init(H.cQpow + i);
            fmpz_pow_ui(c, H.cQ, L);
            _fball_set_fmpz(H.cQpow, c);
            for (i = 1; i < J; i++)
                fball_mul(H.cQpow + i, H.cQpow + i - 1, H.cQpow + i - 1, wp);
            /* the same rule for R: separate its content only when that
               dominates (otherwise R stays short against Q and T) */
            {
                double rc = (double) fmpz_bits(H.cR) - 0.5, rbits;
                rbits = _log2_bound_fmpz(R, Rlen, N) - rc;
                H.have_cR = !fmpz_is_one(H.cR) && rc >= 8.0 * rbits;
            }
            if (H.have_cR)
            {
                fmpz_pow_ui(c, H.cR, L);
                _fball_set_fmpz(H.cRpow, c);
                for (i = 1; i < J; i++)
                    fball_mul(H.cRpow + i, H.cRpow + i - 1, H.cRpow + i - 1, wp);
            }
            fmpz_clear(c);

            bsplit_pow2(T, Qt, Rt, 0, K, J, 0, &H);

            /* Q = Q' cQ^N, cQ^N = (cQ^(N/2))^2 */
            fball_mul(t, H.cQpow + (J - 1), H.cQpow + (J - 1), wp);
            fball_mul(num, Qt, t, wp);
            fball_swap(Qt, num);
            for (i = 0; i < 2 * J; i++)
                fball_clear(H.cQpow + i);
            flint_free(H.cQpow);
        }
        else
        {
            H.have_cR = 0;
            bsplit_range(T, Qt, Rt, 1, N + 1, J, 0, &H);
        }

        for (i = 0; i < 5 * depth; i++)
            fball_clear(H.tmp + i);
        flint_free(H.tmp);
    }

    flint_free(H.buf);

    /* num = coefQ Q + coefP T (+/- |coefP| 2^tau |Q|), den = coefD Q */
    _fball_mul_fmpz(num, T, coefP, wp);
    if (tau > -HUGE_VAL)
    {
        /* |Q| < B^exp */
        double e = tau + (double) fmpz_bits(coefP) + 1.0;
        _fball_add_error_2exp(num, (slong) ceil(e) + FLINT_BITS * Qt->exp);
    }
    if (!fmpz_is_zero(coefQ))
    {
        _fball_mul_fmpz(T, Qt, coefQ, wp);
        fball_add(num, num, T, wp);
    }

    if (fmpz_is_one(coefD))
    {
        fball_swap(den, Qt);
    }
    else
    {
        _fball_mul_fmpz(den, Qt, coefD, wp);
    }

    if (s->power == 1)
        fball_div(res, num, den, wp);
    else
        fball_div(res, den, num, wp);

    fball_clear(T);
    fball_clear(Qt);
    fball_clear(Rt);
    fball_clear(num);
    fball_clear(den);
    fball_clear(t);

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
_int_set_int64(fixed_hypgeom_int_struct * x, nn_ptr d, int64_t v)
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
fball_hypgeom_series_int64(fball_t res, int power, int64_t coefP,
    int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen,
    const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n)
{
    fixed_hypgeom_series_struct s;
    fixed_hypgeom_int_struct * c;
    nn_ptr d;
    slong i, tot = Plen + Qlen + Rlen + 3;

    c = flint_malloc(tot * sizeof(fixed_hypgeom_int_struct));
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

    fball_hypgeom_series(res, &s, n);

    flint_free(c);
    flint_free(d);
}
