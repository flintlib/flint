/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* atan(p/q) and atanh(p/q), 0 < p < q given as mpn integers, by binary
   splitting in mp_real arithmetic (the replacement for
   arb_atan_frac_bsplit in the Machin-type precomputations).

   With x = p/q, P = p^2, Q = q^2, d_k = 2k + 1 and s = -1 for atan
   (s = 1 for atanh),

       atan(x), atanh(x) = x sum_{k>=0} s^k (P/Q)^k / d_k.

   Over a range [a, b) of the series, with D = prod d_k,

       sum_{k=a}^{b-1} s^(k-a) (P/Q)^(k-a) / d_k = N / (D Q^(b-a-1)),

   and ranges compose as

       N = N1 D2 Q^(b-m) + s^(m-a) P^(m-a) N2 D1,    D = D1 D2.

   POWER TABLE.  The factors Q^(b-m) and P^(m-a) are pure powers.  The
   tree is made perfectly balanced -- Nt = K L terms, K a power of
   two, L terms per leaf -- so that the only powers needed are
   Q^(L 2^i) (and P^(L 2^i) when p > 1), one per level, built by
   squarings.  A merge then costs one full-size product (by the power)
   and two short-by-long ones by D, where carrying the full
   denominator DEN = D Q^(b-a) through the tree (N = N1 DEN2 +
   s^(m-a) P^(m-a) N2 D1, DEN = DEN1 DEN2) costs two full-size
   products and one short-by-long.  But D = prod d_k grows by
   lg(2 Nt) bits per term against 2 lg q for Q, so the products by D
   are only short for large q, and the table costs a final Q^Nt: it
   pays off (by 8-12% at 10^6 bits for q >= 2^40) when Q has at least
   four times the bits of the d_k and n >= 64, and the carried
   denominator is used otherwise (3-7% faster for small q or small n).

   LEAVES.  Blocks of L terms are summed exactly by nonallocating mpn
   code in one scratch buffer per call, top down,

       N <- DEN + s P d_j N,    DEN <- d_j Q DEN,    D <- d_j D,

   with DEN = D Q^(b-1-j) local to the leaf: a copy and one fused
   mpn_addmul_1 / mpn_submul_1 per term whenever P d_j fits a limb
   (in particular for p = 1), and one mpn_mul_1 for d_j Q whenever that
   fits (q < 2^32 roughly), multi-limb products otherwise.

   ERRORS.  Merges and the final division are mp_real operations; the
   series tail, at most x^(2 Nt + 1) / ((2 Nt + 1)(1 - x^2)), is added
   as a radius. */

#define ATF_LEAF_LIMBS_MAX 64
#define ATF_LEAF_TERMS_MAX 32

typedef struct
{
    int alt;                /* atan: s = -1 */
    nn_srcptr P, Q;         /* p^2, q^2 */
    slong Pn, Qn;
    int P_one;
    slong L;                /* terms per leaf */
    slong wp;
    nn_ptr buf;             /* leaf scratch: 5 arrays of bufn limbs */
    slong bufn;
    mp_real_struct * Qpow;    /* Qpow[i] = Q^(L 2^i) */
    mp_real_struct * Ppow;    /* Ppow[i] = P^(L 2^i) (unless p = 1) */
    mp_real_struct * tmp;     /* merge temporaries, 5 per level */
    double pt;              /* bits per term of N, DEN */
}
atf_t;

#define MUL1(x, xn, c)                                          \
    do {                                                        \
        ulong __cy = mpn_mul_1((x), (x), *(xn), (c));           \
        (x)[*(xn)] = __cy;                                      \
        *(xn) += (__cy != 0);                                   \
    } while (0)

static slong
_normalize(nn_srcptr x, slong xn)
{
    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    return xn;
}

/* exact N, D (and DEN = D Q^(b-a)) over [a, b) */
static void
leaf(nn_ptr * Np, slong * nlp, nn_ptr * Dp, slong * dlp, nn_ptr * Ep,
    slong * elp, slong a, slong b, const atf_t * A)
{
    slong bufn = A->bufn;
    nn_ptr Nb = A->buf, Wb = Nb + bufn, Db = Wb + bufn, Eb = Db + bufn,
        Tb = Eb + bufn;
    slong nl = 0, wl, dl = 1, el = 1, tl, j;
    ulong cy, dq = 0;
    int dq_fits;

    Db[0] = 1;
    Eb[0] = 1;

    for (j = b - 1; j >= a; j--)
    {
        ulong dj = 2 * (ulong) j + 1;

        /* N <- DEN + s P d_j N, in W as a copy of DEN plus one fused
           pass (P d_j N < DEN: the added part is shorter) */
        flint_mpn_copyi(Wb, Eb, el);
        wl = el;
        if (nl != 0)
        {
            nn_srcptr src = Nb;
            slong sl = nl;
            ulong m = dj;

            if (!A->P_one)
            {
                ulong hi, lo;
                if (A->Pn == 1)
                    umul_ppmm(hi, lo, A->P[0], dj);
                else
                    hi = 1;

                if (hi == 0)
                {
                    m = lo;
                }
                else
                {
                    /* T = N P, then the fused pass by d_j */
                    if (nl >= A->Pn)
                        flint_mpn_mul(Tb, Nb, nl, A->P, A->Pn);
                    else
                        flint_mpn_mul(Tb, A->P, A->Pn, Nb, nl);
                    src = Tb;
                    sl = _normalize(Tb, nl + A->Pn);
                }
            }

            FLINT_ASSERT(wl >= sl);
            if (A->alt)
            {
                cy = mpn_submul_1(Wb, src, sl, m);
                if (wl > sl)
                    cy = mpn_sub_1(Wb + sl, Wb + sl, wl - sl, cy);
                FLINT_ASSERT(cy == 0);
                wl = _normalize(Wb, wl);
            }
            else
            {
                cy = mpn_addmul_1(Wb, src, sl, m);
                if (wl > sl)
                    cy = mpn_add_1(Wb + sl, Wb + sl, wl - sl, cy);
                Wb[wl] = cy;
                wl += (cy != 0);
            }
        }
        FLINT_SWAP(nn_ptr, Nb, Wb);
        nl = wl;

        /* DEN <- d_j Q DEN */
        dq_fits = 0;
        if (A->Qn == 1)
        {
            ulong hi;
            umul_ppmm(hi, dq, A->Q[0], dj);
            dq_fits = (hi == 0);
        }
        if (dq_fits)
        {
            MUL1(Eb, &el, dq);
        }
        else
        {
            if (A->Qn == 1)
            {
                MUL1(Eb, &el, A->Q[0]);
            }
            else if (A->Qn <= 4)
            {
                /* short Q: mpn_mul_1 and mpn_addmul_1 passes, without
                   the dispatch of a general product */
                slong t;
                Tb[el] = mpn_mul_1(Tb, Eb, el, A->Q[0]);
                for (t = 1; t < A->Qn; t++)
                    Tb[el + t] = mpn_addmul_1(Tb + t, Eb, el, A->Q[t]);
                tl = _normalize(Tb, el + A->Qn);
                flint_mpn_copyi(Eb, Tb, tl);
                el = tl;
            }
            else
            {
                if (el >= A->Qn)
                    flint_mpn_mul(Tb, Eb, el, A->Q, A->Qn);
                else
                    flint_mpn_mul(Tb, A->Q, A->Qn, Eb, el);
                tl = _normalize(Tb, el + A->Qn);
                flint_mpn_copyi(Eb, Tb, tl);
                el = tl;
            }
            MUL1(Eb, &el, dj);
        }

        MUL1(Db, &dl, dj);

        FLINT_ASSERT(nl + A->Pn + 2 <= bufn && el + A->Qn + 2 <= bufn);
    }

    *Np = Nb;
    *nlp = nl;
    *Dp = Db;
    *dlp = dl;
    *Ep = Eb;
    *elp = el;
}

static void bsplit_den(mp_real_t N, mp_real_t D, mp_real_t E, slong a0, slong k,
    slong lev, const atf_t * A, int need_d);
static void bsplit(mp_real_t N, mp_real_t D, slong a0, slong k, slong lev,
    const atf_t * A);

/* a subtree as a job for _mp_real_parallel_pair (E == NULL: the power
   table variant); with own set it runs on private leaf buffers and
   level temporaries, freed at the end (the powers are shared, read
   only) */
typedef struct
{
    mp_real_struct * N, * D, * E;
    slong a0, k, lev;
    int need_d, own;
    const atf_t * A;
}
atf_job;

static void
_atf_job(void * arg)
{
    atf_job * J = (atf_job *) arg;
    atf_t A2;
    const atf_t * A = J->A;
    slong i, nt = 5 * (J->lev + 1);

    if (J->own)
    {
        A2 = *J->A;
        A2.buf = flint_malloc(5 * A2.bufn * sizeof(ulong));
        A2.tmp = flint_malloc(nt * sizeof(mp_real_struct));
        for (i = 0; i < nt; i++)
            mp_real_init(A2.tmp + i);
        A = &A2;
    }

    if (J->E != NULL)
        bsplit_den(J->N, J->D, J->E, J->a0, J->k, J->lev, A, J->need_d);
    else
        bsplit(J->N, J->D, J->a0, J->k, J->lev, A);

    if (J->own)
    {
        for (i = 0; i < nt; i++)
            mp_real_clear(A2.tmp + i);
        flint_free(A2.tmp);
        flint_free(A2.buf);
    }
}

/* subtrees of at least ATF_PAR_LEAVES leaves may run on two threads,
   by the rule of MP_REAL_PAR_CAP: 0 serial, 1 the halves on two
   threads, 2 above the cap (the merge on two threads) */
#define ATF_PAR_LEAVES 64

static int
_atf_par(slong k, const atf_t * A)
{
    if (k < ATF_PAR_LEAVES || flint_get_num_threads() < 2)
        return 0;
    return ((double) (k * A->L) * A->pt
        > MP_REAL_PAR_CAP * FLINT_BITS * (double) A->wp) ? 2 : 1;
}

/* the halves of the 2h leaves from a0 at level lev */
static void
_atf_halves(mp_real_struct * N, mp_real_struct * D, mp_real_struct * E,
    mp_real_struct * N2, mp_real_struct * D2, mp_real_struct * E2, slong a0,
    slong h, slong lev, int need_d, int den, int fork, const atf_t * A)
{
    if (fork)
    {
        atf_job L, R;
        L.N = N; L.D = D; L.E = E; L.a0 = a0; L.k = h; L.lev = lev - 1;
        L.need_d = 1; L.own = 0; L.A = A;
        R.N = N2; R.D = D2; R.E = E2; R.a0 = a0 + h; R.k = h;
        R.lev = lev - 1; R.need_d = need_d; R.own = 1; R.A = A;
        _mp_real_parallel_pair(_atf_job, &L, _atf_job, &R);
    }
    else if (den)
    {
        bsplit_den(N, D, E, a0, h, lev - 1, A, 1);
        bsplit_den(N2, D2, E2, a0 + h, h, lev - 1, A, need_d);
    }
    else
    {
        bsplit(N, D, a0, h, lev - 1, A);
        bsplit(N2, D2, a0 + h, h, lev - 1, A);
    }
}

/* the merges as two independent halves (on two threads above the cap):
   with the carried denominator {T = N1 DEN2, DEN = DEN1 DEN2} and
   {U = N2 D1 P^h, D = D1 D2}, with the power table {T = N1 D2 Q^(L h)}
   and {U = N2 D1 P^h, D = D1 D2}; then N = T +/- U */
typedef struct
{
    mp_real_struct * N, * D, * E, * N2, * D2, * E2, * T, * U;
    const mp_real_struct * qpow, * ppow;
    int need_d;
    slong wp;
}
atf_merge_struct;

static void
_atf_merge_x(void * arg)
{
    atf_merge_struct * M = (atf_merge_struct *) arg;
    if (M->E != NULL)
    {
        mp_real_mul(M->T, M->N, M->E2, M->wp);
        mp_real_mul(M->E, M->E, M->E2, M->wp);
    }
    else
    {
        mp_real_mul(M->T, M->N, M->D2, M->wp);
        mp_real_mul(M->T, M->T, M->qpow, M->wp);
    }
}

static void
_atf_merge_y(void * arg)
{
    atf_merge_struct * M = (atf_merge_struct *) arg;
    mp_real_mul(M->U, M->N2, M->D, M->wp);
    if (M->ppow != NULL)
        mp_real_mul(M->U, M->U, M->ppow, M->wp);
    if (M->need_d)
        mp_real_mul(M->D, M->D, M->D2, M->wp);
}

static void
_atf_merge(atf_merge_struct * M, int sub, int par)
{
    if (par)
        _mp_real_parallel_pair(_atf_merge_x, M, _atf_merge_y, M);
    else
    {
        _atf_merge_x(M);
        _atf_merge_y(M);
    }
    if (sub)
        mp_real_sub(M->N, M->T, M->U, M->wp);
    else
        mp_real_add(M->N, M->T, M->U, M->wp);
}

/* free this level's temporaries when large: they are otherwise held
   through all the merges above */
#define ATF_KEEP 2048

static void
_atf_release(mp_real_struct * T)
{
    slong i;
    if (T[0].alloc > ATF_KEEP || T[1].alloc > ATF_KEEP)
        for (i = 0; i < 5; i++)
        {
            mp_real_clear(T + i);
            mp_real_init(T + i);
        }
}

/* the variant carrying the full denominator DEN = D Q^(len) through
   the tree: N = N1 DEN2 + s^h P^h N2 D1, DEN = DEN1 DEN2, D = D1 D2
   (not needed along the right spine); temporaries: 5 per level */
static void
bsplit_den(mp_real_t N, mp_real_t D, mp_real_t E, slong a0, slong k, slong lev,
    const atf_t * A, int need_d)
{
    if (k == 1)
    {
        nn_ptr Nb, Db, Eb;
        slong nl, dl, el;

        leaf(&Nb, &nl, &Db, &dl, &Eb, &el, a0 * A->L, (a0 + 1) * A->L, A);
        _mp_real_set_mpn_2exp(N, Nb, nl, 0);
        _mp_real_set_mpn_2exp(D, Db, dl, 0);
        _mp_real_set_mpn_2exp(E, Eb, el, 0);
    }
    else
    {
        mp_real_struct * N2 = A->tmp + 5 * lev, * D2 = N2 + 1, * T = N2 + 2,
            * U = N2 + 3, * E2 = N2 + 4;
        slong h = k / 2, wp = A->wp;

        int par = _atf_par(k, A);
        atf_merge_struct M;

        _atf_halves(N, D, E, N2, D2, E2, a0, h, lev, need_d, 1, par == 1, A);

        M.N = N; M.D = D; M.E = E; M.N2 = N2; M.D2 = D2; M.E2 = E2;
        M.T = T; M.U = U; M.qpow = NULL;
        M.ppow = A->P_one ? NULL : A->Ppow + (lev - 1);
        M.need_d = need_d; M.wp = wp;
        _atf_merge(&M, A->alt && ((A->L * h) & 1), par == 2);
        _atf_release(N2);
    }
}

/* the range of k leaves starting at leaf index a0 (k a power of two),
   level lev = log2(k), with the power table; temporaries of level lev
   at A->tmp + 5 lev */
static void
bsplit(mp_real_t N, mp_real_t D, slong a0, slong k, slong lev,
    const atf_t * A)
{
    if (k == 1)
    {
        nn_ptr Nb, Db, Eb;
        slong nl, dl, el;

        leaf(&Nb, &nl, &Db, &dl, &Eb, &el, a0 * A->L, (a0 + 1) * A->L, A);
        _mp_real_set_mpn_2exp(N, Nb, nl, 0);
        _mp_real_set_mpn_2exp(D, Db, dl, 0);
    }
    else
    {
        mp_real_struct * N2 = A->tmp + 5 * lev, * D2 = N2 + 1, * T = N2 + 2,
            * U = N2 + 3;
        slong h = k / 2, wp = A->wp;
        mp_real_struct * qpow = A->Qpow + (lev - 1);   /* Q^(L h) */

        int par = _atf_par(k, A);
        atf_merge_struct M;

        _atf_halves(N, D, NULL, N2, D2, NULL, a0, h, lev, 1, 0, par == 1, A);

        M.N = N; M.D = D; M.E = NULL; M.N2 = N2; M.D2 = D2; M.E2 = NULL;
        M.T = T; M.U = U; M.qpow = qpow;
        M.ppow = A->P_one ? NULL : A->Ppow + (lev - 1);
        M.need_d = 1; M.wp = wp;
        _atf_merge(&M, A->alt && ((A->L * h) & 1), par == 2);
        _atf_release(N2);
    }
}

/* pow[i] = base^(L 2^i), i < J: base^L exactly, then squarings */
static void
_powtab(mp_real_struct * pow, nn_srcptr base, slong bn, slong L, slong J,
    slong n)
{
    nn_ptr w, v, w0;
    slong wl, vl, t, i;

    w0 = w = flint_malloc(2 * (L * bn + 2) * sizeof(ulong));
    v = w + (L * bn + 2);

    flint_mpn_copyi(w, base, bn);
    wl = bn;
    for (t = 1; t < L; t++)
    {
        if (wl >= bn)
            flint_mpn_mul(v, w, wl, base, bn);
        else
            flint_mpn_mul(v, base, bn, w, wl);
        vl = _normalize(v, wl + bn);
        FLINT_SWAP(nn_ptr, v, w);
        wl = vl;
    }
    _mp_real_set_mpn_2exp(pow, w, wl, 0);
    for (i = 1; i < J; i++)
        mp_real_mul(pow + i, pow + i - 1, pow + i - 1, n);

    flint_free(w0);
}

/* lower and upper bounds for log2 of (x, xn), xn >= 1, top limb
   nonzero */
static void
_log2_bounds(double * lo, double * hi, nn_srcptr x, slong xn)
{
    double t = (double) x[xn - 1];
    slong e = FLINT_BITS * (xn - 1);

    if (xn >= 2)
    {
        t = t * ldexp(1.0, FLINT_BITS) + (double) x[xn - 2];
        e -= FLINT_BITS;
    }

    /* t is within a relative 2^-52 of the top two limbs, and x / 2^e
       lies below them plus one unit when lower limbs were dropped
       (the log2 of doubles is accurate to a few ulps, covered by the
       2^-50 margins) */
    *lo = log2(t * (1.0 - 0x1p-50)) + (double) e;
    *hi = log2((t + (xn > 2 ? 1.0 : 0.0)) * (1.0 + 0x1p-50)) + (double) e;
}

static void
_atan_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q,
    slong qn, int hyperbolic, slong n)
{
    atf_t A;
    double lp_lo, lp_hi, lq_lo, lq_hi, lx, pt, tail;
    slong Nt0, Nt, K, J, L, Lmax, i, tb;
    nn_ptr P, Q, pq;
    slong Pn, Qn, pqn;
    mp_real_t N, D, T;
    TMP_INIT;

    pn = _normalize(p, pn);
    qn = _normalize(q, qn);

    /* two guard limbs of our own: a limb count is only a precision up
       to the size of the top limb, and a deep tree compounds errors */
    n += 2;

    if (pn == 0)
    {
        mp_real_zero(res);
        return;
    }

    if (qn < pn || (qn == pn && mpn_cmp(p, q, pn) >= 0))
        flint_throw(FLINT_ERROR, "_mp_real_atan_frac_bsplit: need p < q\n");

    TMP_START;

    _log2_bounds(&lp_lo, &lp_hi, p, pn);
    _log2_bounds(&lq_lo, &lq_hi, q, qn);
    /* lower bound for lg(q/p); the series is only meant for p/q
       bounded away from 1 (x <= 0.99 here) */
    lx = lq_lo - lp_hi;
    if (!(lx >= 0.015))
        flint_throw(FLINT_ERROR,
            "_mp_real_atan_frac_bsplit: p/q too close to 1\n");

    /* terms for a tail x^(2 Nt + 1) / (1 - x^2) below 2^-tb */
    tb = FLINT_BITS * n + 16;
    Nt0 = (slong) ceil((tb / lx - 1.0) / 2.0) + 1;
    Nt0 = FLINT_MAX(Nt0, 1);

    P = TMP_ALLOC((2 * pn + 2 * qn + pn + qn) * sizeof(ulong));
    Q = P + 2 * pn;
    pq = Q + 2 * qn;
    flint_mpn_sqr(P, p, pn);
    Pn = _normalize(P, 2 * pn);
    flint_mpn_sqr(Q, q, qn);
    Qn = _normalize(Q, 2 * qn);
    if (pn >= qn)
        flint_mpn_mul(pq, p, pn, q, qn);
    else
        flint_mpn_mul(pq, q, qn, p, pn);
    pqn = _normalize(pq, pn + qn);

    A.alt = !hyperbolic;
    A.P = P;
    A.Pn = Pn;
    A.P_one = (Pn == 1 && P[0] == 1);
    A.Q = Q;
    A.Qn = Qn;
    A.wp = n;

    /* leaf size: about min(wp, 64) limbs of numerator bits, at most 32
       terms; then K = 2^J leaves of L terms, Nt = K L >= Nt0 */
    pt = FLINT_BITS * (Qn - 1) + FLINT_BIT_COUNT(Q[Qn - 1])
        + FLINT_BIT_COUNT(2 * (ulong) Nt0 + 1) + 1;
    Lmax = (slong) ((FLINT_BITS * FLINT_MIN(n, ATF_LEAF_LIMBS_MAX)) / pt);
    Lmax = FLINT_MAX(Lmax, 1);
    Lmax = FLINT_MIN(Lmax, ATF_LEAF_TERMS_MAX);
    for (J = 0, K = 1; K * Lmax < Nt0; J++, K *= 2)
        ;
    L = (Nt0 + K - 1) / K;
    Nt = K * L;
    A.L = L;

    A.bufn = (slong) ((L + 2) * pt) / FLINT_BITS + Pn + Qn + 8;
    A.pt = pt;

    mp_real_init(N);
    mp_real_init(D);
    mp_real_init(T);

    A.buf = flint_malloc(5 * A.bufn * sizeof(ulong));

    if (K == 1)
    {
        /* a single exact leaf: x N / (D Q^(Nt-1)) = p q N / DEN */
        nn_ptr Nb, Db, Eb;
        slong nl, dl, el;

        leaf(&Nb, &nl, &Db, &dl, &Eb, &el, 0, Nt, &A);
        _mp_real_set_mpn_2exp(N, Nb, nl, 0);
        _mp_real_set_mpn_2exp(D, Eb, el, 0);
        _mp_real_set_mpn_2exp(T, pq, pqn, 0);
        mp_real_mul(N, N, T, n);
        mp_real_div(res, N, D, n);
    }
    else
    {
        /* the power table pays off when Q is long against the d_k,
           i.e. when the products by D are short (measured: 8-12%
           faster for q >= 2^40 at 10^6 bits, 3-7% slower for small q
           or at low precision, where D is comparable to Q^len and the
           final Q^Nt costs more than the merge product it saves) */
        int powtab = (FLINT_BITS * (Qn - 1) + FLINT_BIT_COUNT(Q[Qn - 1])
            >= 4 * FLINT_BIT_COUNT(2 * (ulong) Nt + 1)) && n >= 64;

        A.Qpow = flint_malloc(2 * J * sizeof(mp_real_struct));
        A.Ppow = A.Qpow + J;
        A.tmp = flint_malloc(5 * (J + 1) * sizeof(mp_real_struct));
        for (i = 0; i < 2 * J; i++)
            mp_real_init(A.Qpow + i);
        for (i = 0; i < 5 * (J + 1); i++)
            mp_real_init(A.tmp + i);

        if (!A.P_one)
            _powtab(A.Ppow, P, Pn, L, J, n);

        if (powtab)
        {
            _powtab(A.Qpow, Q, Qn, L, J, n);
            bsplit(N, D, 0, K, J, &A);

            /* p q N / (D Q^Nt), Q^Nt = (Q^(Nt/2))^2 */
            mp_real_mul(T, A.Qpow + (J - 1), A.Qpow + (J - 1), n);
            mp_real_mul(D, D, T, n);
        }
        else
        {
            /* p q N / DEN */
            mp_real_t E;
            mp_real_init(E);
            bsplit_den(N, D, E, 0, K, J, &A, 0);
            mp_real_swap(D, E);
            mp_real_clear(E);
        }

        _mp_real_set_mpn_2exp(T, pq, pqn, 0);
        mp_real_mul(N, N, T, n);
        mp_real_div(res, N, D, n);

        for (i = 0; i < 2 * J; i++)
            mp_real_clear(A.Qpow + i);
        for (i = 0; i < 5 * (J + 1); i++)
            mp_real_clear(A.tmp + i);
        flint_free(A.Qpow);
        flint_free(A.tmp);
    }

    flint_free(A.buf);

    /* tail: x^(2 Nt + 1) / ((2 Nt + 1)(1 - x^2)) <= 2^tail */
    {
        double xu = exp2(-lx);      /* upper bound for x < 0.99 */
        tail = -(2.0 * Nt + 1.0) * lx - log2(1.0 - xu * xu);
        tail = ceil(tail * (1.0 + 1e-12) + 1e-9) + 1;
    }
    {
        slong t = (slong) tail;
        slong qq = t >> (FLINT_BITS == 64 ? 6 : 5);
        _mp_real_add_error_ulps_at(res, ldexp(1.0, (int) (t - qq * FLINT_BITS)), qq);
    }

    mp_real_clear(N);
    mp_real_clear(D);
    mp_real_clear(T);
    TMP_END;
}

void
_mp_real_atan_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q,
    slong qn, slong n)
{
    _atan_frac_bsplit(res, p, pn, q, qn, 0, n);
}

void
_mp_real_atanh_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q,
    slong qn, slong n)
{
    _atan_frac_bsplit(res, p, pn, q, qn, 1, n);
}
