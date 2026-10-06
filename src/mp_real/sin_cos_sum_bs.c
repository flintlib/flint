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

/* Joint binary splitting for the sine and cosine series of
   x' = x B^-D (B = 2^64, D a LIMB count) over their COMMON
   denominator: with y = x^2 and N terms of the y-series,

       sin(x')     = x' (Q B^Qexp - B_ B^be) / (Q B^Qexp),
       1 - cos(x') = A B^ae / (Q B^Qexp),

   where Q B^Qexp collects the products (2k)(2k+1) together with the
   per-term frame B^(2Dk) (whole zero limbs of Q stripped into
   Qexp at every merge) and

       -B_ / (Q B^Qexp) = sum_{k=1}^{N} (-y)^k
                                        / (prod (2i)(2i+1) B^(2Dk)),
       -A / (Q B^Qexp) = sum_{k=1}^{N} (2k+1) (-y)^k
                                        / (prod (2i)(2i+1) B^(2Dk))

   -- the second series is 1 - cos rewritten over (2k+1)! via
   1/(2k)! = (2k+1)/(2k+1)!, which is what makes ONE tree with ONE
   denominator serve both components: the per-term weights (2k+1)
   are absolute, so they compose across merges unchanged, and the
   A-track rides along for two extra multiplications per merge
   (5 in all against 3 for a single series) -- half the work of a
   Gaussian splitting of exp(ix) over all 2N terms, which this
   replaced (measured 1.3x slower than arb's per-slice sine
   splitting; see dev/notes).

   Signs: every subrange sum T(a,b) = sum_{k=a+1}^{b} (-1)^(k-a) |t_k|
   has a NEGATIVE leading term, so both tracks are stored as
   magnitudes U = -T; the merge

       U(a,b) = U(a,m) Q2 B^(Q2exp) + (-1)^(step+1) y^step U(m,b)

   subtracts when step = m - a is odd and adds when it is even,
   dominance-safe because the subtrahend carries y^step B^(-2D step)
   times a ratio of index products below 1 (asserted).  The
   bridging denominator Q(a,m) B^(QE(a,m)) absorbs the whole
   B^(-2D step) exactly, whatever the limb split between mantissa
   and exponent, so no explicit frame appears in the merge.

   The power table (over y) and the split-step exponents mirror
   exp_sum_bs.c / the arb original; instead of the bit version's
   per-leaf 2-adic folding, whole zero LIMBS of Q are stripped into
   Qexp at every merge -- the same content to within one limb, at
   limb granularity, and the table entries strip their dead low
   limbs the same way (sparse arguments in particular).

   TRUNCATION: the exact splitting integers span about
   (b - a)(2 x 64 D + eps) bits -- several times the target
   precision at the root -- but the burst factor keeps only a
   cap-limb window, so nothing below relative 2^(-64 lmax) can
   matter.  Both numerator
   tracks (and the power-table entries) therefore carry an explicit
   LSB exponent and are truncated to lmax limbs whenever they grow
   past it, each drop one-sided (understating) below one ulp; the
   giant top-of-tree multiplications then run at lmax limbs instead
   of the full exact span, the former full-span shifts become
   exponent arithmetic, and the per-merge scratch stops scaling
   with the subtree.  Q stays exact: its bit length
   ~ N lg N < 64 lmax throughout the burst's range.  Over <= lg N
   truncated levels the accumulated relative error stays far below
   the factor window's own 2^(-64 cap) truncation, inside the same
   budget.  All frames (ae, be, Qexp, the table exponents) are LIMB
   counts, so every alignment inside a merge is a limb move -- the
   finest frame that keeps the result within lmax + 1 limbs of its
   top is chosen, dropping under one ulp of the result -- and no
   bit shift occurs anywhere in the kernel. */

/* exponents (m - a) reachable from [0, n) with m = a + (b-a)/2 and
   b - a = 2 inlined; duplicated from exp_sum_bs.c (static there) */
static slong
compute_bs_exponents(slong * tab, slong n)
{
    slong a, b, aa, ba, bb, length;

    if (n == 1)
    {
        tab[0] = 1;
        return 1;
    }
    if (n == 2 || n == 3 || n == 4)
    {
        tab[0] = 1;
        tab[1] = 2;
        return 2;
    }
    if (n == 6)
    {
        tab[0] = 1;
        tab[1] = 2;
        tab[2] = 3;
        return 3;
    }

    a = n >> 1;
    b = n - (n >> 1);
    tab[0] = a;
    length = 1;

    for (;;)
    {
        aa = a >> 1;
        ba = b >> 1;
        bb = b - ba;

        tab[length] = ba;
        length++;

        if (ba == 3)
        {
            tab[length] = 2;
            tab[length + 1] = 1;
            length += 2;
            break;
        }

        if (ba == 1 || (ba == 2 && (n & (n - 1)) == 0))
            break;

        if (aa != ba && aa != 1)
        {
            tab[length] = aa;
            length++;
        }

        a = aa;
        b = bb;
    }

    if (tab[length - 1] != 1)
    {
        tab[length] = 1;
        length++;
    }

    for (a = 0; a < length / 2; a++)
    {
        b = tab[a];
        tab[a] = tab[length - a - 1];
        tab[length - a - 1] = b;
    }

    return length;
}

static slong
get_exp_pos(const slong * tab, slong step)
{
    slong i;
    for (i = 0; ; i++)
    {
        if (tab[i] == step)
            return i;
        FLINT_ASSERT(tab[i] != 0);
    }
}

/* ---- small mpn-number helpers: (ptr, len), len 0 means zero ---- */

static slong
nnn_normalize(nn_srcptr p, slong n)
{
    while (n > 0 && p[n - 1] == 0)
        n--;
    return n;
}

static slong
nnn_mul(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    if (an == 0 || bn == 0)
        return 0;
    if (an >= bn)
        flint_mpn_mul(z, a, an, b, bn);
    else
        flint_mpn_mul(z, b, bn, a, an);
    return nnn_normalize(z, an + bn);
}

static slong
nnn_add(nn_ptr z, slong zn, nn_srcptr a, slong an)
{
    ulong cy;
    if (an == 0)
        return zn;
    if (zn >= an)
    {
        cy = mpn_add(z, z, zn, a, an);
        z[zn] = cy;
        return zn + (cy != 0);
    }
    flint_mpn_zero(z + zn, an - zn);
    cy = mpn_add(z, a, an, z, zn);
    z[an] = cy;
    return an + (cy != 0);
}

/* z -= a in place, dominance-guaranteed z > a; returns new length */
static slong
nnn_sub(nn_ptr z, slong zn, nn_srcptr a, slong an)
{
    if (an == 0)
        return zn;
    FLINT_ASSERT(zn >= an);
    {
        ulong bw = mpn_sub(z, z, zn, a, an);
        FLINT_ASSERT(bw == 0);
        (void) bw;
    }
    return nnn_normalize(z, zn);
}

/* strip trailing zero LIMBS of q into the limb exponent *e (the
   denominator's unstripped even factors carry ~2N bits of 2-adic
   content across a range of N terms, whole dead limbs at scale) */
static slong
nnn_strip_low(nn_ptr q, slong qn, slong * e)
{
    slong z = 0;
    while (z < qn - 1 && q[z] == 0)
        z++;
    if (z > 0)
    {
        flint_mpn_copyi(q, q + z, qn - z);
        *e += z;
        qn -= z;
    }
    return qn;
}

/* ---- the recursion ---- */

typedef struct
{
    const slong * xexp;
    nn_srcptr * xpow;
    const slong * xlen;
    const slong * xe;       /* LSB LIMB exponents of the table */
    slong D;                /* argument frame in limbs: x B^-D,
                               the y = x^2 series scaled B^(-2D)
                               per term */
    slong lmax;             /* mantissa truncation length (limbs) */
}
bs_args;

static slong
qbound(const bs_args * args, slong a, slong b)
{
    return ((b - a) * 2 * FLINT_BIT_COUNT(2 * (ulong) b + 1))
        / FLINT_BITS + 3;
}

/* room for an A or B track over [a, b): the exact value spans at most
   (b - a) terms of 2D limbs plus its factorial content, and the
   truncation caps it at lmax (the scratch of a merge was sized by
   lmax alone at first, which made every node of the tree, down to
   the leaves, malloc and touch a full-precision buffer: at 3 10^7 bits
   a third of the tree's time went to those page faults) */
static slong
ebound(const bs_args * args, slong a, slong b)
{
    return FLINT_MIN(args->lmax, (b - a) * 2 * args->D + 2)
        + qbound(args, a, b) + 4;
}

/* normalize (z, zn) B^(ze): first strip trailing zero LIMBS into
   the exponent (exact -- the limb-granular 2-valuation rule, which
   also catches sparse arguments whose powers carry dead low
   limbs), then truncate to at most args->lmax limbs by dropping
   low limbs into the exponent, one-sided below one ulp */
static slong
nnn_trunc(const bs_args * args, nn_ptr z, slong zn, slong * ze)
{
    slong d = 0;
    while (d < zn - 1 && z[d] == 0)
        d++;
    if (d > 0)
    {
        flint_mpn_copyi(z, z + d, zn - d);
        *ze += d;
        zn -= d;
    }
    if (zn > args->lmax)
    {
        d = zn - args->lmax;
        flint_mpn_copyi(z, z + d, args->lmax);
        *ze += d;
        return args->lmax;
    }
    return zn;
}

/* z = z B^(ze) +- p B^(pe): both operands are aligned to the
   FINEST limb frame that keeps the result within lmax + 1 limbs of
   its top -- the exact frame min(ze, pe) whenever that fits, and
   otherwise raised just enough, so that any dropped limbs sit at
   least lmax limbs below the result's top (a one-sided relative
   2^(-64 lmax + 63), the same class as the explicit truncations).
   All frames are limb counts, so every alignment is a limb move --
   no bit shifts.  sub = 1 subtracts (dominance-guaranteed, so no
   cancellation of the top).  z must have room for lmax + 4 limbs
   plus the incoming length; tmp for lmax + 5.  Returns the new
   length; *ze gets the result frame. */
static slong
nnn_combine(const bs_args * args, nn_ptr z, slong zn, slong * ze,
    nn_srcptr p, slong pn, slong pe, int sub, nn_ptr tmp)
{
    slong tb, f, q, pl;

    if (pn == 0)
        return zn;
    if (zn == 0)
    {
        FLINT_ASSERT(!sub);
        flint_mpn_copyi(z, p, pn);
        *ze = pe;
        return pn;
    }

    /* top-limb bound of the result and the target frame */
    tb = FLINT_MAX(*ze + zn, pe + pn) + 1;
    f = FLINT_MIN(*ze, pe);
    if (f < tb - (args->lmax + 1))
        f = tb - (args->lmax + 1);

    /* align z to frame f in place */
    if (*ze > f)
    {
        q = *ze - f;
        memmove(z + q, z, zn * sizeof(ulong));
        flint_mpn_zero(z, q);
        zn = zn + q;
    }
    else if (*ze < f)
    {
        q = f - *ze;
        if (q >= zn)
        {
            z[0] = 0;
            zn = 0;
        }
        else
        {
            flint_mpn_copyi(z, z + q, zn - q);
            zn = zn - q;
        }
    }
    *ze = f;

    /* align p to frame f in tmp (or subtract in place when its
       frame sits at or above f: offset addition/subtraction) */
    if (pe >= f)
    {
        q = pe - f;
        flint_mpn_zero(tmp, q);
        flint_mpn_copyi(tmp + q, p, pn);
        pl = pn + q;
    }
    else
    {
        q = f - pe;
        if (q >= pn)
            pl = 0;
        else
        {
            flint_mpn_copyi(tmp, p + q, pn - q);
            pl = nnn_normalize(tmp, pn - q);
        }
    }

    if (pl == 0)
        return FLINT_MAX(zn, 1);
    if (zn == 0)
    {
        FLINT_ASSERT(!sub);
        flint_mpn_copyi(z, tmp, pl);
        return pl;
    }
    if (sub)
        zn = nnn_sub(z, zn, tmp, pl);
    else
        zn = nnn_add(z, zn, tmp, pl);
    return zn;
}

/* (t, returned length) = (x, xn) (2a + 4)(2a + 5): one word when the
   product fits (always on 64-bit machines), else two multiplications;
   t has room for xn + 2 limbs */
static slong
_mul_f1(nn_ptr t, nn_srcptr x, slong xn, slong a)
{
    ulong hi, f, a1 = (ulong) (2 * a + 4), a2 = (ulong) (2 * a + 5);
    slong l;

    umul_ppmm(hi, f, a1, a2);
    if (hi == 0)
    {
        t[xn] = mpn_mul_1(t, x, xn, f);
        return nnn_normalize(t, xn + 1);
    }
    t[xn] = mpn_mul_1(t, x, xn, a1);
    l = nnn_normalize(t, xn + 1);
    t[l] = mpn_mul_1(t, t, l, a2);
    return nnn_normalize(t, l + 1);
}

static void bsplit(nn_ptr A, slong * an, slong * ae, nn_ptr B, slong * bn,
    slong * be, nn_ptr Q, slong * qn, slong * QE,
    const bs_args * args, slong a, slong b);

/* THREADS.  A subtree of at least TRIG_BS_PAR_TERMS terms whose
   merge products reach TRIG_BS_PAR_LIMBS limbs runs its two halves on
   two threads (_mp_real_parallel_pair, which splits the budget between
   them for the forks below), and then its merge as two jobs, the A
   track and the B track with the Q product, each on its own scratch.
   Below the thresholds, or without threads, the serial code runs. */
#ifndef TRIG_BS_PAR_LIMBS
#define TRIG_BS_PAR_LIMBS 2048
#endif
#ifndef TRIG_BS_PAR_TERMS
#define TRIG_BS_PAR_TERMS 16
#endif

typedef struct
{
    nn_ptr A, B, Q;
    slong * an, * ae, * bn, * be, * qn, * QE;
    const bs_args * args;
    slong a, b;
}
trig_bs_job;

static void
_trig_bs_job(void * arg)
{
    trig_bs_job * J = (trig_bs_job *) arg;
    bsplit(J->A, J->an, J->ae, J->B, J->bn, J->be, J->Q, J->qn, J->QE,
        J->args, J->a, J->b);
}

/* the merge as two jobs: with the A track, {A Q2, P A2} and
   {B Q2, P B2, Q Q2}; without it (the sine-only trees), {P B2} and
   {B Q2, Q Q2} -- the two products of a track being independent, the
   combine follows the jobs */
typedef struct
{
    nn_ptr A, B, Q, sca, scpa, scb, scpb, scq;
    slong * an, * ae, * bn, * be, * qn, * QE;
    nn_srcptr A2, B2, Q2, P;
    slong a2n, b2n, q2n, Q2exp, Pl;
    slong la, lpa, lb, lpb, lq;
    const bs_args * args;
}
trig_merge_struct;

static void
_trig_merge_x(void * arg)
{
    trig_merge_struct * M = (trig_merge_struct *) arg;
    if (M->A != NULL)
    {
        M->la = nnn_mul(M->sca, M->A, *M->an, M->Q2, M->q2n);
        M->lpa = nnn_mul(M->scpa, M->P, M->Pl, M->A2, M->a2n);
    }
    else
        M->lpb = nnn_mul(M->scpb, M->P, M->Pl, M->B2, M->b2n);
}

static void
_trig_merge_y(void * arg)
{
    trig_merge_struct * M = (trig_merge_struct *) arg;
    M->lb = nnn_mul(M->scb, M->B, *M->bn, M->Q2, M->q2n);
    if (M->A != NULL)
        M->lpb = nnn_mul(M->scpb, M->P, M->Pl, M->B2, M->b2n);
    M->lq = nnn_mul(M->scq, M->Q, *M->qn, M->Q2, M->q2n);
}

static void
bsplit(nn_ptr A, slong * an, slong * ae, nn_ptr B, slong * bn,
    slong * be, nn_ptr Q, slong * qn, slong * QE,
    const bs_args * args, slong a, slong b)
{
    if (b - a == 1)
    {
        {
            ulong hi;
            umul_ppmm(hi, Q[0], (ulong) (2 * a + 2),
                (ulong) (2 * a + 3));
            Q[1] = hi;
            *qn = 1 + (hi != 0);
        }
        *QE = 2 * args->D;

        flint_mpn_copyi(B, args->xpow[0], args->xlen[0]);
        *bn = args->xlen[0];
        *be = args->xe[0];
        *bn = nnn_trunc(args, B, *bn, be);
        if (A != NULL)
        {
            A[args->xlen[0]] = mpn_mul_1(A, args->xpow[0],
                args->xlen[0], (ulong) (2 * a + 3));
            *an = nnn_normalize(A, args->xlen[0] + 1);
            *ae = args->xe[0];
            *an = nnn_trunc(args, A, *an, ae);
        }
    }
    else if (b - a == 2)
    {
        /* U(a, a+2) = t1 B^(2D) - t2 in each track, all at the
           y-table frame; the B^(2D) offset goes into a frame gap
           handled by nnn_combine as a pure limb move */
        nn_ptr t1, t2, tt;
        slong l1, l2, e1;
        TMP_INIT;
        TMP_START;
        /* the leaf values are a few limbs (y and y^2 with their
           factors), and the frame gap 2D of the combine */
        t1 = TMP_ALLOC((args->xlen[0] + args->xlen[1] * 2 + 2 * args->D + 24)
            * sizeof(ulong));
        t2 = t1 + args->xlen[0] + 6;
        tt = t2 + args->xlen[1] + 6;

        /* B */
        l1 = _mul_f1(t1, args->xpow[0], args->xlen[0], a);
        e1 = args->xe[0] + 2 * args->D;
        flint_mpn_copyi(B, t1, l1);
        *bn = l1;
        *be = e1;
        *bn = nnn_combine(args, B, *bn, be, args->xpow[1],
            args->xlen[1], args->xe[1], 1, tt);
        *bn = nnn_trunc(args, B, *bn, be);

        /* A */
        if (A != NULL)
        {
            l1 = _mul_f1(t1, args->xpow[0], args->xlen[0], a);
            t1[l1] = mpn_mul_1(t1, t1, l1, (ulong) (2 * a + 3));
            l1 = nnn_normalize(t1, l1 + 1);
            flint_mpn_copyi(A, t1, l1);
            *an = l1;
            *ae = e1;
            t2[args->xlen[1]] = mpn_mul_1(t2, args->xpow[1],
                args->xlen[1], (ulong) (2 * a + 5));
            l2 = nnn_normalize(t2, args->xlen[1] + 1);
            *an = nnn_combine(args, A, *an, ae, t2, l2, args->xe[1],
                1, tt);
            *an = nnn_trunc(args, A, *an, ae);
        }
        TMP_END;

        {
            ulong q1[2], q2[2], hi;
            slong ln1, ln2;
            umul_ppmm(hi, q1[0], (ulong) (2 * a + 2),
                (ulong) (2 * a + 3));
            q1[1] = hi;
            ln1 = 1 + (hi != 0);
            umul_ppmm(hi, q2[0], (ulong) (2 * a + 4),
                (ulong) (2 * a + 5));
            q2[1] = hi;
            ln2 = 1 + (hi != 0);
            *QE = 4 * args->D;
            *qn = nnn_mul(Q, q1, ln1, q2, ln2);
        }
    }
    else if (b - a >= TRIG_BS_PAR_TERMS && flint_get_num_threads() >= 2
        && args->lmax >= TRIG_BS_PAR_LIMBS)
    {
        slong step, m, i, a2n, b2n, q2n, a2e, b2e, pe, Q2exp;
        slong scl, ttl, e2l;
        nn_ptr A2, B2, Q2, sca, scpa, scb, scpb, scq, tt;
        trig_bs_job L, R;
        trig_merge_struct M;
        TMP_INIT;

        step = (b - a) / 2;
        m = a + step;

        /* the five products of the merge each in their own scratch:
           X Q2 up to lmax + q2n limbs, P X2 up to 2 lmax, Q Q2 within
           qbound(a, b); the right half's tracks sized for [m, b) */
        scl = 2 * ebound(args, a, b) + 8;
        ttl = ebound(args, a, b) + 6;
        e2l = ebound(args, m, b) + qbound(args, m, b) + 4;

        TMP_START;
        A2 = TMP_ALLOC((2 * e2l + qbound(args, m, b) + 5 * scl + ttl)
            * sizeof(ulong));
        B2 = A2 + e2l;
        Q2 = B2 + e2l;
        sca = Q2 + qbound(args, m, b);
        scpa = sca + scl;
        scb = scpa + scl;
        scpb = scb + scl;
        scq = scpb + scl;
        tt = scq + scl;

        L.A = A; L.an = an; L.ae = ae; L.B = B; L.bn = bn; L.be = be;
        L.Q = Q; L.qn = qn; L.QE = QE; L.args = args; L.a = a; L.b = m;
        R.A = (A == NULL) ? NULL : A2; R.an = &a2n; R.ae = &a2e;
        R.B = B2; R.bn = &b2n; R.be = &b2e; R.Q = Q2; R.qn = &q2n;
        R.QE = &Q2exp; R.args = args; R.a = m; R.b = b;
        _mp_real_parallel_pair(_trig_bs_job, &L, _trig_bs_job, &R);

        i = get_exp_pos(args->xexp, step);
        pe = args->xe[i];

        M.A = A; M.an = an; M.ae = ae; M.B = B; M.bn = bn; M.be = be;
        M.Q = Q; M.qn = qn; M.QE = QE;
        M.sca = sca; M.scpa = scpa; M.scb = scb; M.scpb = scpb;
        M.scq = scq;
        M.A2 = A2; M.a2n = a2n; M.B2 = B2; M.b2n = b2n;
        M.Q2 = Q2; M.q2n = q2n; M.Q2exp = Q2exp;
        M.P = args->xpow[i]; M.Pl = args->xlen[i];
        M.args = args;
        _mp_real_parallel_pair(_trig_merge_x, &M, _trig_merge_y, &M);

        /* A = (A Q2) B^Q2exp -+ P A2, B likewise, truncated to lmax */
        if (A != NULL)
        {
            flint_mpn_copyi(A, sca, M.la);
            *an = M.la;
            *ae += Q2exp;
            *an = nnn_combine(args, A, *an, ae, scpa, M.lpa, pe + a2e,
                (int) (step & 1), tt);
            *an = nnn_trunc(args, A, *an, ae);
        }
        flint_mpn_copyi(B, scb, M.lb);
        *bn = M.lb;
        *be += Q2exp;
        *bn = nnn_combine(args, B, *bn, be, scpb, M.lpb, pe + b2e,
            (int) (step & 1), tt);
        *bn = nnn_trunc(args, B, *bn, be);

        /* Q = Q Q2, exact */
        flint_mpn_copyi(Q, scq, M.lq);
        *qn = M.lq;
        *QE = *QE + Q2exp;
        *qn = nnn_strip_low(Q, *qn, QE);
        if (mp_real_get_verbose())
            _mp_real_progress("trig bsplit: merged terms [%wd, %wd), B has %wd limbs", a, b, *bn);
        TMP_END;
    }
    else
    {
        slong step, m, i, a2n, b2n, q2n, l, a2e, b2e, pe;
        slong Q2exp;
        nn_ptr A2, B2, Q2, sc, tt;
        nn_srcptr P;
        slong Pl;
        TMP_INIT;

        step = (b - a) / 2;
        m = a + step;

        TMP_START;
        /* sc holds one product at a time: A1 Q2 (lmax + qbound) or
           P A2 (both factors up to lmax: 2 lmax limbs).  The A2/B2
           slots are sized for the UNTRUNCATED intermediate
           A Q2 / B Q2 (lmax + q2n limbs) that transiently lives in
           the caller's corresponding slot one level up -- sizing
           them at lmax + 4 alone lets that copy spill into the Q
           slot and clobber the left child's accumulated Q before
           the Q update reads it */
        {
            slong e2l = ebound(args, m, b) + qbound(args, m, b) + 4;
            slong scl = 2 * ebound(args, a, b) + 8;
            A2 = TMP_ALLOC((2 * e2l + qbound(args, m, b) + scl
                + ebound(args, a, b) + 6) * sizeof(ulong));
            B2 = A2 + e2l;
            Q2 = B2 + e2l;
            sc = Q2 + qbound(args, m, b);
            tt = sc + scl;
        }

        bsplit(A, an, ae, B, bn, be, Q, qn, QE, args, a, m);
        bsplit((A == NULL) ? NULL : A2, &a2n, &a2e, B2, &b2n,
            &b2e, Q2, &q2n, &Q2exp, args, m, b);

        i = get_exp_pos(args->xexp, step);
        P = args->xpow[i];
        Pl = args->xlen[i];
        pe = args->xe[i];

        /* A = (A Q2) 2^(Q2exp) -+ y^step A2, frames aligned by
           nnn_combine, truncated to lmax; skipped in the sine-only
           mode (A == NULL) */
        if (A != NULL)
        {
            l = nnn_mul(sc, A, *an, Q2, q2n);
            flint_mpn_copyi(A, sc, l);
            *an = l;
            *ae += Q2exp;
            l = nnn_mul(sc, P, Pl, A2, a2n);
            *an = nnn_combine(args, A, *an, ae, sc, l, pe + a2e,
                (int) (step & 1), tt);
            *an = nnn_trunc(args, A, *an, ae);
        }

        l = nnn_mul(sc, B, *bn, Q2, q2n);
        flint_mpn_copyi(B, sc, l);
        *bn = l;
        *be += Q2exp;
        l = nnn_mul(sc, P, Pl, B2, b2n);
        *bn = nnn_combine(args, B, *bn, be, sc, l, pe + b2e,
            (int) (step & 1), tt);
        *bn = nnn_trunc(args, B, *bn, be);

        /* Q = Q Q2, exact */
        l = nnn_mul(sc, Q, *qn, Q2, q2n);
        flint_mpn_copyi(Q, sc, l);
        *qn = l;
        *QE = *QE + Q2exp;
        *qn = nnn_strip_low(Q, *qn, QE);
        if (*bn >= MP_REAL_VERBOSE_MIN_LIMBS && mp_real_get_verbose())
            _mp_real_progress("trig bsplit: merged terms [%wd, %wd), B has %wd limbs", a, b, *bn);
        TMP_END;
    }
}

/* A (length *an) B^(*ae), B (*bn) B^(*be), Q (*qn), QE as
   documented in the file comment, all exponents LIMB counts: over
   N terms of the y = x^2 series,
   1 - cos(x B^-D) ~ A B^ae / (Q B^QE) and
   sin(x B^-D) ~ (x B^-D) (Q B^QE - B B^be) / (Q B^QE),
   with A and B truncated to lmax limbs (one-sided low; the
   working accuracy class is ~2^45 ulps at lmax across the tree's
   drops, absorbed by the caller's guard limbs).  A and B must each
   have room for (lmax + (N 2 bits(2N+1)) / 64 + 8) limbs and Q for
   ((N 2 bits(2N+1)) / 64 + 3) limbs.  Requires xn <= D (value
   below 1) and lmax >= 4. */
static void
_sum_bs_powtab_core(nn_ptr A, slong * an, slong * ae,
    nn_ptr B, slong * bn, slong * be, nn_ptr Q, slong * qn,
    slong * QE, nn_srcptr xp, slong xn,
    slong D, slong N, slong lmax)
{
    slong xexp[2 * FLINT_BITS];
    slong xlen[2 * FLINT_BITS];
    slong xe[2 * FLINT_BITS];
    nn_srcptr xpow[2 * FLINT_BITS];
    nn_ptr storage, y;
    slong yn, ye, length, i, off;
    bs_args args;

    FLINT_ASSERT(N >= 1 && xn >= 1 && xp[xn - 1] != 0);
    FLINT_ASSERT(xn <= D);      /* value x B^-D below 1 */
    /* dominance in the merge subtractions needs y B^-2D well below
       one: D >= 1 keeps a huge margin for any realistic N */
    FLINT_ASSERT(D >= 1);
    FLINT_ASSERT(lmax >= 4);

    for (i = 0; i < 2 * FLINT_BITS; i++)
        xexp[i] = 0;
    length = compute_bs_exponents(xexp, N);

    args.lmax = lmax;
    args.D = D;

    /* y = x^2; the power table is over y, entries truncated to
       lmax limbs with LSB exponents (each truncation one-sided
       below one ulp; y^k accumulates under k of them, absorbed by
       the same budget as the merge truncations) */
    y = flint_malloc(2 * xn * sizeof(ulong));
    flint_mpn_sqr(y, xp, xn);
    yn = nnn_normalize(y, 2 * xn);
    ye = 0;
    yn = nnn_trunc(&args, y, yn, &ye);

    {
        slong total = 0;
        for (i = 0; i < length; i++)
            total += FLINT_MIN(xexp[i] * 2 * xn + 1, lmax + 2);
        storage = flint_malloc(total * sizeof(ulong));
    }
    off = 0;
    {
        /* every entry is built in scratch (a square of a truncated
           entry spans up to 2 lmax limbs) and stored truncated */
        nn_ptr sq = flint_malloc((4 * lmax + 8) * sizeof(ulong));
        for (i = 0; i < length; i++)
        {
            nn_ptr dst = storage + off;
            slong sl, se, src;
            off += FLINT_MIN(xexp[i] * 2 * xn + 1, lmax + 2);

            if (i == 0)
            {
                FLINT_ASSERT(xexp[0] == 1);
                flint_mpn_copyi(dst, y, yn);
                xlen[0] = yn;
                xe[0] = ye;
                xpow[0] = dst;
                continue;
            }
            if (xexp[i] == 2 * xexp[i - 1])
                src = i - 1;
            else if (i >= 2 && xexp[i] == 2 * xexp[i - 2])
                src = i - 2;
            else if (xexp[i] == 2 * xexp[i - 1] + 1)
                src = i - 1;
            else
            {
                FLINT_ASSERT(i >= 2
                    && xexp[i] == 2 * xexp[i - 2] + 1);
                src = i - 2;
            }

            flint_mpn_sqr(sq, xpow[src], xlen[src]);
            sl = nnn_normalize(sq, 2 * xlen[src]);
            se = 2 * xe[src];
            sl = nnn_trunc(&args, sq, sl, &se);

            if (xexp[i] & 1)
            {
                /* odd: multiply the square by y */
                sl = nnn_mul(sq + lmax + 2, sq, sl, y, yn);
                se += ye;
                flint_mpn_copyi(sq, sq + lmax + 2, sl);
                sl = nnn_trunc(&args, sq, sl, &se);
            }
            flint_mpn_copyi(dst, sq, sl);
            xlen[i] = sl;
            xe[i] = se;
            xpow[i] = dst;
        }
        flint_free(sq);
    }

    args.xexp = xexp;
    args.xpow = xpow;
    args.xlen = xlen;
    args.xe = xe;

    bsplit(A, an, ae, B, bn, be, Q, qn, QE, &args, 0, N);

    flint_free(storage);
    flint_free(y);
}

void
_mp_real_sin_cos_sum_bs_powtab(nn_ptr A, slong * an, slong * ae,
    nn_ptr B, slong * bn, slong * be, nn_ptr Q, slong * qn,
    slong * QE, nn_srcptr xp, slong xn,
    slong D, slong N, slong lmax)
{
    _sum_bs_powtab_core(A, an, ae, B, bn, be, Q, qn, QE, xp, xn,
        D, N, lmax);
}
