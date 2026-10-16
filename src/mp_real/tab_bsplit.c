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

/* Binary splitting for the table values log(1 + 2^-i) and
   atan(2^-i) in mp_real arithmetic.

   Both are sums V = sum_{k>=0} s^k |r|^k / d_k, d_k = dbase + dstep k,
   with either an alternating dyadic ratio -2^-e or (for
   log(1 + 2^-i) = 2 atanh(1/q), q = 2^(i+1) + 1) the ratio q^-2:

       atan(2^-i)      = 2^-i  sum (-1)^k 2^(-2ik) / (2k+1)
       log(1 + 2^-i)   = 2^-i  sum (-1)^k 2^(-ik)  / (k+1)
                       = (2/q) sum q^(-2k) / (2k+1).

   Over a range [a, b), let U(a, b) = sum_{k=a}^{b-1} s^(k-a)
   |r|^(k-a+1) / d_k, D(a, b) = prod d_k and DEN(a, b) = D |r|^-(b-a),
   and write U = N / DEN.  Ranges then compose with

       N = N1 DEN2 + s^(m-a) N2 D1,   DEN = DEN1 DEN2,   D = D1 D2,

   and V = U(0, Nt) / |r|.  For a dyadic ratio DEN2 = D2 2^(e(b-m)) is
   a shift and never stored (three products per merge); for q^-2 it
   is carried (four products per merge, D skipped along the right
   spine where no left sibling needs it).

   LEAVES.  Blocks of L terms are evaluated exactly by nonallocating
   mpn code in one scratch buffer per call.  For the dyadic ratio, N
   is summed directly as N = sum (-1)^(k-a) (D / d_k) 2^(e (b-1-k)):
   pieces as short as D = prod d_k (a few limbs) at their own shifts,
   so a leaf costs O(L |D|), not O(L |N|).  For q^-2 the Horner scheme

       N <- DEN(j+1) + d_j N,   D <- d_j D,   DEN <- d_j q^2 DEN

   costs a copy, one fused mpn_addmul_1 and single-limb
   multiplications per term (d_j q^2 is one mpn_mul_1 while it fits a
   limb).  A leaf spans about min(working precision, 64 limbs) of
   numerator bits and at most 32 terms.  When the whole series fits in
   one leaf of at most 64 terms and twice that size (small precision
   or large i), no mp_real arithmetic is needed at all: the exact
   partial sum is floored by one division.

   ERRORS.  Everything above the leaves is mp_real arithmetic, which
   tracks the truncation errors (and the series tail, added as a
   radius) rigorously; the output is a lower bound read off the ball,
   retrying at higher precision if the radius were ever too wide.
   Either way the result is floor(V B^n) or one below it (the tail
   is < 2^-24 ulps and the alternating series use an even term count,
   so the exact partial sums never exceed V). */

/* tuning parameters: leaf size bounds, the exact single-leaf range
   (terms, and size relative to a leaf), the direct-series cutoff for
   the logarithm (i >= C1, or i >= C2 when n >= C3), guard limbs */
#define TAB_BS_LEAF_LIMBS_MAX 64
#define TAB_BS_LEAF_TERMS_MAX 32
#define TAB_BS_EXACT_TERMS_MAX 64
#define TAB_BS_EXACT_FACTOR 2
#define TAB_BS_LOG_C1 30
#define TAB_BS_LOG_C2 20
#define TAB_BS_LOG_C3 32
#define TAB_BS_GUARD 3

typedef struct
{
    int atanh;          /* ratio q^-2 (else -2^-e) */
    slong e;            /* dyadic shift per term */
    ulong i;
    ulong dbase, dstep; /* d_k = dbase + dstep k */
    slong L;            /* terms per leaf */
    slong wp;           /* working precision, limbs */
    nn_ptr buf;         /* leaf scratch: 6 arrays of bufn limbs */
    slong bufn;
    mp_real_struct * tmp; /* merge temporaries: 4 per tree level, reused
                           by all nodes of a level (never live at the
                           same time), so allocation happens once */
}
tab_bs_t;

/* (x, *xn) *= c */
#define MUL1(x, xn, c)                                          \
    do {                                                        \
        ulong __cy = mpn_mul_1((x), (x), *(xn), (c));           \
        (x)[*(xn)] = __cy;                                      \
        *(xn) += (__cy != 0);                                   \
    } while (0)

/* (x, *xn) += (y, yn); x must have room for max + 1 limbs */
static void
_add(nn_ptr x, slong * xn, nn_srcptr y, slong yn)
{
    ulong cy;

    if (*xn >= yn)
    {
        cy = mpn_add(x, x, *xn, y, yn);
    }
    else
    {
        cy = mpn_add(x, y, yn, x, *xn);
        *xn = yn;
    }
    x[*xn] = cy;
    *xn += (cy != 0);
}

/* (t, *tn) = (x, xn) << sh, xn >= 1 with top limb nonzero */
static void
_shl(nn_ptr t, slong * tn, nn_srcptr x, slong xn, slong sh)
{
    slong q = sh / FLINT_BITS;
    int r = sh % FLINT_BITS;

    flint_mpn_zero(t, q);
    if (r)
    {
        ulong cy = mpn_lshift(t + q, x, xn, r);
        t[q + xn] = cy;
        *tn = q + xn + (cy != 0);
    }
    else
    {
        flint_mpn_copyi(t + q, x, xn);
        *tn = q + xn;
    }
}

/* scratch layout: six arrays of bufn limbs */
#define BUF_N0(P) ((P)->buf)
#define BUF_N1(P) ((P)->buf + (P)->bufn)
#define BUF_D(P)  ((P)->buf + 2 * (P)->bufn)
#define BUF_E(P)  ((P)->buf + 3 * (P)->bufn)
#define BUF_T(P)  ((P)->buf + 4 * (P)->bufn)
#define BUF_S(P)  ((P)->buf + 5 * (P)->bufn)

/* dyadic leaf as a sum of short pieces: with D = prod d_k,

       N = sum_{k=a}^{b-1} (-1)^(k-a) (D / d_k) 2^(e (b-1-k)),

   each piece as long as D (a few limbs: the d_k are small) at its own
   shift, so the cost is O((b - a) |D|) rather than the O((b - a) |N|)
   of a Horner scheme, whose every step touches all of N -- the
   difference is large when e is.  Even and odd k accumulate
   separately and are subtracted once. */
static void
leaf_dyadic(nn_ptr * Np, slong * nlp, slong * dlp, slong a,
    slong b, const tab_bs_t * P)
{
    nn_ptr Pb = BUF_N0(P), Mb = BUF_N1(P), Db = BUF_D(P), Qb = BUF_T(P),
        Sb = BUF_S(P);
    slong dl = 1, ql, len, k, q;
    ulong cy;

    Db[0] = 1;
    for (k = a; k < b; k++)
        MUL1(Db, &dl, P->dbase + P->dstep * (ulong) k);

    /* the top piece D / d_a sits at shift e (b-1-a) */
    len = (P->e * (b - 1 - a)) / FLINT_BITS + dl + 2;
    FLINT_ASSERT(len + 1 <= P->bufn);
    flint_mpn_zero(Pb, len);
    flint_mpn_zero(Mb, len);

    for (k = a; k < b; k++)
    {
        slong sh = P->e * (b - 1 - k);
        nn_ptr acc = ((k - a) & 1) ? Mb : Pb;

        mpn_divexact_1(Qb, Db, dl, P->dbase + P->dstep * (ulong) k);
        ql = dl;
        while (Qb[ql - 1] == 0)
            ql--;

        q = sh / FLINT_BITS;
        if (sh % FLINT_BITS)
        {
            Sb[ql] = mpn_lshift(Sb, Qb, ql, sh % FLINT_BITS);
            cy = mpn_add(acc + q, acc + q, len - q, Sb, ql + 1);
        }
        else
        {
            cy = mpn_add(acc + q, acc + q, len - q, Qb, ql);
        }
        FLINT_ASSERT(cy == 0);
        (void) cy;
    }

    /* positive: the terms decrease */
    cy = mpn_sub_n(Pb, Pb, Mb, len);
    FLINT_ASSERT(cy == 0);
    while (Pb[len - 1] == 0)
        len--;

    *Np = Pb;
    *nlp = len;
    *dlp = dl;
}

/* atanh leaf by the Horner scheme */
static void
leaf_atanh(nn_ptr * Np, slong * nlp, slong * dlp, slong * elp, slong a,
    slong b, const tab_bs_t * P, int need_d)
{
    nn_ptr Nb = BUF_N0(P), Wb = BUF_N1(P), Db = BUF_D(P), Eb = BUF_E(P),
        Tb = BUF_T(P), Sb = BUF_S(P);
    slong nl = 0, wl, dl = 1, el = 1, tl, sl, j;
    ulong cy;

    Db[0] = 1;
    Eb[0] = 1;

    for (j = b - 1; j >= a; j--)
    {
        ulong dj = P->dbase + P->dstep * (ulong) j;

        /* N <- DEN + d_j N, formed in W as a copy of DEN plus one
           fused pass (N < DEN: the added part is shorter) */
        flint_mpn_copyi(Wb, Eb, el);
        wl = el;
        if (nl != 0)
        {
            FLINT_ASSERT(wl >= nl);
            cy = mpn_addmul_1(Wb, Nb, nl, dj);
            if (wl > nl)
                cy = mpn_add_1(Wb + nl, Wb + nl, wl - nl, cy);
            Wb[wl] = cy;
            wl += (cy != 0);
        }
        FLINT_SWAP(nn_ptr, Nb, Wb);
        nl = wl;

        /* DEN <- d_j q^2 DEN, one mpn_mul_1 when d_j q^2 fits */
        if (2 * P->i + 2 < FLINT_BITS)
        {
            ulong q2 = (UWORD(1) << (2 * P->i + 2))
                + (UWORD(1) << (P->i + 2)) + 1, hi, lo;

            umul_ppmm(hi, lo, q2, dj);
            if (hi == 0)
            {
                MUL1(Eb, &el, lo);
            }
            else
            {
                MUL1(Eb, &el, q2);
                MUL1(Eb, &el, dj);
            }
        }
        else
        {
            _shl(Tb, &tl, Eb, el, 2 * P->i + 2);
            _shl(Sb, &sl, Eb, el, P->i + 2);
            _add(Tb, &tl, Sb, sl);
            _add(Tb, &tl, Eb, el);
            flint_mpn_copyi(Eb, Tb, tl);
            el = tl;
            MUL1(Eb, &el, dj);
        }

        if (need_d)
            MUL1(Db, &dl, dj);

        FLINT_ASSERT(nl + 2 <= P->bufn && el + 2 <= P->bufn
            && dl + 2 <= P->bufn);
    }

    *Np = Nb;
    *nlp = nl;
    *dlp = dl;
    *elp = el;
}

/* exact N, D, DEN over [a, b): N at *Np (one of BUF_N0, BUF_N1),
   D at BUF_D, DEN at BUF_E (atanh only) */
static void
leaf_mpn(nn_ptr * Np, slong * nlp, slong * dlp, slong * elp, slong a,
    slong b, const tab_bs_t * P, int need_d)
{
    if (P->atanh)
    {
        leaf_atanh(Np, nlp, dlp, elp, a, b, P, need_d);
    }
    else
    {
        leaf_dyadic(Np, nlp, dlp, a, b, P);
        *elp = 0;
    }
}

static void
leaf(mp_real_t N, mp_real_t D, mp_real_t E, slong a, slong b,
    const tab_bs_t * P, int need_d)
{
    nn_ptr Nb;
    slong nl, dl, el;

    leaf_mpn(&Nb, &nl, &dl, &el, a, b, P, need_d);

    _mp_real_set_mpn_2exp(N, Nb, nl, 0);
    if (need_d || !P->atanh)
        _mp_real_set_mpn_2exp(D, BUF_D(P), dl, 0);
    if (P->atanh)
        _mp_real_set_mpn_2exp(E, BUF_E(P), el, 0);
}

static void
bsplit(mp_real_t N, mp_real_t D, mp_real_t E, slong a, slong b,
    const tab_bs_t * P, int need_d, slong depth)
{
    mp_real_struct * N2 = P->tmp + 4 * depth, * D2 = N2 + 1, * E2 = N2 + 2,
        * T = N2 + 3;
    slong m, wp = P->wp;

    if (b - a <= P->L)
    {
        leaf(N, D, E, a, b, P, need_d);
        return;
    }

    m = a + P->L * (((b - a + P->L - 1) / P->L) / 2);

    bsplit(N, D, E, a, m, P, 1, depth + 1);
    bsplit(N2, D2, E2, m, b, P, need_d || !P->atanh, depth + 1);

    /* outputs never alias inputs (mp_real would stage the product in
       a temporary and copy it back); results are swapped into place */
    if (!P->atanh)
    {
        /* N = N1 D2 2^(e (b-m)) + s^(m-a) N2 D1, D = D1 D2 */
        mp_real_mul(T, N, D2, wp);
        mp_real_mul_2exp_si(T, T, P->e * (b - m));
        mp_real_mul(E2, N2, D, wp);
        if ((m - a) & 1)
            mp_real_sub(N, T, E2, wp);
        else
            mp_real_add(N, T, E2, wp);
        mp_real_mul(T, D, D2, wp);
        mp_real_swap(D, T);
    }
    else
    {
        /* N = N1 DEN2 + N2 D1, DEN = DEN1 DEN2, D = D1 D2 */
        mp_real_mul(T, N, E2, wp);
        mp_real_mul(N, N2, D, wp);
        mp_real_add(N2, T, N, wp);
        mp_real_swap(N, N2);
        mp_real_mul(T, E, E2, wp);
        mp_real_swap(E, T);
        if (need_d)
        {
            mp_real_mul(T, D, D2, wp);
            mp_real_swap(D, T);
        }
    }
}

/* res = floor(X 2^s / Y), known to be < B^n */
static void
_div_floor(nn_ptr res, slong n, nn_srcptr X, slong xn, nn_srcptr Y,
    slong yn, slong s)
{
    nn_ptr num, den, q, r;
    slong num_n, den_n, qn;
    TMP_INIT;

    TMP_START;

    if (s >= 0)
    {
        num = TMP_ALLOC((xn + s / FLINT_BITS + 1) * sizeof(ulong));
        _shl(num, &num_n, X, xn, s);
        den = (nn_ptr) Y;
        den_n = yn;
    }
    else
    {
        den = TMP_ALLOC((yn - s / FLINT_BITS + 1) * sizeof(ulong));
        _shl(den, &den_n, Y, yn, -s);
        num = (nn_ptr) X;
        num_n = xn;
    }

    flint_mpn_zero(res, n);

    if (num_n >= den_n)
    {
        qn = num_n - den_n + 1;
        q = TMP_ALLOC((qn + den_n) * sizeof(ulong));
        r = q + qn;
        flint_mpn_tdiv_qr(q, r, num, num_n, den, den_n);
        while (qn > 0 && q[qn - 1] == 0)
            qn--;
        FLINT_ASSERT(qn <= n);
        flint_mpn_copyi(res, q, qn);
    }

    TMP_END;
}

/* enclose also an additive perturbation of magnitude <= 2^t */
static void
_mp_real_add_error_2exp(mp_real_t x, slong t)
{
    slong q = t >> (FLINT_BITS == 64 ? 6 : 5);    /* floor */
    _mp_real_add_error_ulps_at(x, (double) (UWORD(1) << (t - q * FLINT_BITS)), q);
}

/* res = n-limb lower bound for v in [0, 1): at most the true value
   and, when the radius is below B/16 guard ulps, at most one ulp
   below its floor; returns 0 if the radius is wider */
static int
_get_lower(nn_ptr res, slong n, const mp_real_t v)
{
    nn_ptr y;
    ulong err;
    int ok = 0;
    TMP_INIT;

    TMP_START;
    y = TMP_ALLOC((n + 1) * sizeof(ulong));
    _mp_real_get_fixed(y, &err, v, n + 1);

    if (err < (UWORD(1) << (FLINT_BITS - 4)))
    {
        /* |true - y| <= err guard ulps */
        if (mpn_sub_1(y, y, n + 1, err + 1))
            flint_mpn_zero(y, n + 1);
        flint_mpn_copyi(res, y + 1, n);
        ok = 1;
    }

    TMP_END;
    return ok;
}

static void
_tab_bs_setup(tab_bs_t * P, slong guard, slong n, slong Nt, slong pt)
{
    slong leafbits;

    P->wp = n + guard;
    leafbits = FLINT_BITS * FLINT_MIN(P->wp, TAB_BS_LEAF_LIMBS_MAX);
    P->L = FLINT_MAX(leafbits / pt, 1);
    P->L = FLINT_MIN(P->L, TAB_BS_LEAF_TERMS_MAX);
    P->L = FLINT_MIN(P->L, Nt);
    P->bufn = (P->L * pt + pt + 2 * FLINT_BITS) / FLINT_BITS + 4;
}

static void
_tab_bs(nn_ptr res, ulong i, slong n, int atanh, slong e, ulong dbase,
    ulong dstep)
{
    tab_bs_t P;
    slong Nt, tb, pt, guard, depth, k;
    mp_real_t N, D, E, v;

    FLINT_ASSERT(i >= 1);
    FLINT_ASSERT(n >= 1);

    /* tail below 2^-tb = 2^-24 output ulps */
    tb = FLINT_BITS * n + 24;

    if (atanh)
    {
        /* tail of (2/q) sum_{k >= Nt} q^(-2k)/(2k+1)
           <= 2^(2 - (i+1)(2 Nt + 1)) */
        Nt = ((tb + 2) / ((slong) i + 1) + 1) / 2 + 1;
    }
    else
    {
        /* alternating decreasing: tail <= 2^-i 2^(-e Nt); an even
           Nt makes the tail nonnegative (used by the exact path) */
        Nt = (tb - (slong) i) / e + 1;
        Nt += (Nt & 1);
    }
    Nt = FLINT_MAX(Nt, 1);

    P.atanh = atanh;
    P.e = e;
    P.i = i;
    P.dbase = dbase;
    P.dstep = dstep;

    /* bits of growth per term of the leaf numerators */
    pt = (atanh ? 2 * (slong) i + 2 : e)
        + FLINT_BIT_COUNT(dbase + dstep * (ulong) Nt) + 1;

    if (Nt <= TAB_BS_EXACT_TERMS_MAX && Nt * pt <= TAB_BS_EXACT_FACTOR
            * FLINT_BITS * FLINT_MAX(n + TAB_BS_GUARD, TAB_BS_LEAF_LIMBS_MAX))
    {
        /* one exact leaf: the floor of the partial sum, which is at
           most the true value since the tail is nonnegative, and
           short by less than one ulp plus the tail */
        slong nl, dl, el, tl;
        nn_ptr Nb, Tb;
        TMP_INIT;

        P.L = Nt;
        P.bufn = (Nt * pt + pt + 2 * FLINT_BITS) / FLINT_BITS + 4;

        TMP_START;
        P.buf = TMP_ALLOC(6 * P.bufn * sizeof(ulong));
        leaf_mpn(&Nb, &nl, &dl, &el, 0, Nt, &P, !atanh);

        if (!atanh)
        {
            /* 2^-i N / (D 2^(e (Nt-1))) */
            _div_floor(res, n, Nb, nl, BUF_D(&P), dl,
                FLINT_BITS * n - (slong) i - e * (Nt - 1));
        }
        else
        {
            /* 2 q N / DEN */
            Tb = BUF_T(&P);
            _shl(Tb, &tl, Nb, nl, (slong) i + 1);
            _add(Tb, &tl, Nb, nl);
            _div_floor(res, n, Tb, tl, BUF_E(&P), el,
                FLINT_BITS * n + 1);
        }

        TMP_END;
        return;
    }

    mp_real_init(N);
    mp_real_init(D);
    mp_real_init(E);
    mp_real_init(v);

    for (guard = TAB_BS_GUARD; ; guard *= 2)
    {
        _tab_bs_setup(&P, guard, n, Nt, pt);
        /* depth: the leaf count halves (rounding up) per level */
        for (depth = 1, k = (Nt + P.L - 1) / P.L; k > 1; k = (k + 1) / 2)
            depth++;
        P.buf = flint_malloc(6 * P.bufn * sizeof(ulong));
        P.tmp = flint_malloc(4 * depth * sizeof(mp_real_struct));
        for (k = 0; k < 4 * depth; k++)
            mp_real_init(P.tmp + k);

        bsplit(N, D, E, 0, Nt, &P, 0, 0);

        for (k = 0; k < 4 * depth; k++)
            mp_real_clear(P.tmp + k);
        flint_free(P.tmp);
        flint_free(P.buf);

        if (!atanh)
        {
            /* V = N / (D 2^(e (Nt-1))); value = 2^-i V */
            mp_real_div(v, N, D, P.wp);
            mp_real_mul_2exp_si(v, v, -((slong) i + e * (Nt - 1)));
            _mp_real_add_error_2exp(v, -((slong) i + e * Nt));
        }
        else
        {
            /* V = q^2 N / DEN; value = (2/q) V = 2 q N / DEN */
            if (i + 1 < FLINT_BITS)
            {
                mp_real_mul_ui(N, N, (UWORD(1) << (i + 1)) + 1, P.wp);
            }
            else
            {
                mp_real_set(v, N);
                mp_real_mul_2exp_si(v, v, (slong) i + 1);
                mp_real_add(N, N, v, P.wp);
            }
            mp_real_div(v, N, E, P.wp);
            mp_real_mul_2exp_si(v, v, 1);
            _mp_real_add_error_2exp(v, 2 - ((slong) i + 1) * (2 * Nt + 1));
        }

        if (_get_lower(res, n, v))
            break;
    }

    mp_real_clear(N);
    mp_real_clear(D);
    mp_real_clear(E);
    mp_real_clear(v);
}

void
_mp_real_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
{
    /* d_k = 2k + 1, ratio -2^(-2i) */
    _tab_bs(res, i, n, 0, 2 * (slong) i, 1, 2);
}

void
_mp_real_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
{
    /* the direct series has twice the terms of the atanh one, but
       dyadic merges (three products instead of four) and cheaper
       leaves; measured faster from i = 20 at moderate precision and
       from i = 30 throughout */
    if (i >= TAB_BS_LOG_C1 || (i >= TAB_BS_LOG_C2 && n >= TAB_BS_LOG_C3))
        _tab_bs(res, i, n, 0, (slong) i, 1, 1);  /* d_k = k + 1 */
    else
        _tab_bs(res, i, n, 1, 0, 1, 2);          /* d_k = 2k + 1 */
}
