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
#include "mpn_extras.h"
#include "fixed.h"

/* Chudnovsky:

       1/pi = 12 sum_{k>=0} (-1)^k (6k)! (A0 + A1 k)
                    / ((3k)! (k!)^3 640320^(3k + 3/2)),

   A0 = 13591409, A1 = 545140134.  As a hypergeometric sum,

       S = sum_{k>=0} A(k) prod_{j<=k} p(j)/q(j),
       A(k) = A0 + A1 k,
       p(k) = -(6k-5)(2k-1)(6k-1),   q(k) = k^3 C,   (k >= 1)
       p(0) = q(0) = 1,              C = 640320^3 / 24,

   whence pi = 640320^(3/2) / (12 S) = D (Q/T) / sqrt(640320) with
   D = 640320^2 / 12 and the binary splitting quantities

       P(a,b) = prod_{k=a}^{b-1} p(k),
       Q(a,b) = prod_{k=a}^{b-1} q(k),
       T(a,b) = Q(a,b) sum_{k=a}^{b-1} A(k) prod_{j=a}^{k} p(j)/q(j)

   satisfying T(a,b) = T(a,m) Q(m,b) + P(a,m) T(m,b), and
   S = T(0,N)/Q(0,N) up to the series tail.

   NORMALIZATION BAKED INTO q(0).  Setting q(0) = c instead of 1
   multiplies Q(0,N) by c while leaving T(0,N) invariant: the k = 0
   term of the sum picks up the compensating 1/c, so
   T = Q sum = (c Q') (sum'/c).  With c = D = 640320^2/12 the
   reconstruction collapses to pi = Q rsqrt(640320) / T with no
   final scalar multiplication, and with c = D/4 the same tree
   yields pi/4 directly -- no trailing bit shift either.

   All quantities are fballs truncated to the caller's precision n:
   exact integers (err = 0.0) until they outgrow n limbs, balls
   afterwards, so the top of the tree runs entirely in O(n)-limb
   windowed arithmetic while the accumulated ulp radius (a few per
   node, ~N/BLK in total, far below the 2^69 normalization threshold
   for any feasible N) is absorbed by the guard limbs in n.

   BASECASE.  Blocks of up to PI_BS_BLK terms are accumulated
   iteratively over exact mpn integers, backward in k with
   Qs = Q(k+1,b) and S = S(k+1) = T(k+1,b) as state:

       S(k) = p(k) (A(k) Qs + S(k+1)),      Qs <- q(k) Qs,

   two limb multiplications and a subtraction per term (S(k) < 0 for
   every k >= 1, since the leading term of the alternating tail
   dominates, so the magnitude recursion is |S| <- g(k) (A(k) Qs -
   |S|), g = -p).  The block product G = prod g(k) accumulates with
   one more mul_1. */

#ifndef PI_BS_BLK
#define PI_BS_BLK 32
#endif

#define PI_A0 UWORD(13591409)
#define PI_A1 UWORD(545140134)
/* C = 640320^3 / 24 = 2^15 3^2 5^3 23^3 29^3, D = 640320^2 / 12 =
   2^10 3 5^2 23^2 29^2, in factor pairs below one 32-bit limb so
   that every width multiplies by whole limbs */
#define PI_C1 UWORD(36864000)              /* 2^15 3^2 5^3 */
#define PI_C2 UWORD(296740963)             /* 23^3 29^3 */
#define PI_D1 UWORD(76800)                 /* 2^10 3 5^2 */
#define PI_D2 UWORD(444889)                /* 23^2 29^2 */
#define PI_D4_1 UWORD(19200)               /* 2^8 3 5^2: (D/4)/D2 */

/* bits per term: |p/q| < 72/C = 2^-47.11 */
#define PI_BITS_PER_TERM 47

/* stack sizes for one basecase block: PI_BS_BLK terms of
   |q(k)| = k^3 C < 2^(54 + 3 log2 kmax) bits, plus slack */
#ifndef PI_BS_ALLOC
#define PI_BS_ALLOC ((PI_BS_BLK * 132) / FLINT_BITS + 12)
#endif

static void
_fball_set_mpn(fball_t x, nn_srcptr d, slong len, int negative)
{
    slong t = 0;

    while (len > 0 && d[len - 1] == 0)
        len--;
    while (t < len && d[t] == 0)
        t++;
    fball_fit(x, len - t);
    flint_mpn_copyi(x->d, d + t, len - t);
    x->size = len - t;
    x->exp = len;
    x->negative = negative && (len > t);
    x->err = 0.0;
}

/* exact P, Q, T over [a, b), 1 <= b - a <= PI_BS_BLK; the leftmost
   leaf (a == 0) multiplies q0f (two whole-limb factors, 1 = skip)
   into q(0) */
static void
pi_bsplit_basecase(fball_t P, fball_t Q, fball_t T, slong a, slong b,
    int need_p, const ulong * q0f)
{
    ulong Qs[PI_BS_ALLOC], S[PI_BS_ALLOC], G[PI_BS_ALLOC],
          tmp[PI_BS_ALLOC];
    slong lq, ls, lg, lt, k;

    /* each linear factor of g(k) = (6k-5)(2k-1)(6k-1) and of
       q(k)/C = k^3 must fit a limb; on 64-bit machines g itself
       fits below PI_FUSE_K and is applied in one multiplication */
    FLINT_ASSERT((ulong) b <= (UWORD_MAX - 1) / 6);

    Qs[0] = 1; lq = 1;
    G[0] = 1; lg = 1;
    ls = 0;

    for (k = b - 1; k >= a; k--)
    {
        ulong cy;

        /* tmp = A(k) Qs - |S|  (positive: the A Qs term
           dominates); A(k) = A0 + A1 k needs two limbs once
           A1 k overflows (always on 32-bit machines beyond
           small k) */
        if ((ulong) k <= (UWORD_MAX - PI_A0) / PI_A1)
        {
            ulong A = PI_A0 + PI_A1 * (ulong) k;
            cy = mpn_mul_1(tmp, Qs, lq, A);
            tmp[lq] = cy;
            lt = lq + (cy != 0);
        }
        else
        {
            ulong ahi, alo;
            umul_ppmm(ahi, alo, PI_A1, (ulong) k);
            add_ssaaaa(ahi, alo, ahi, alo, UWORD(0), PI_A0);
            tmp[lq] = mpn_mul_1(tmp, Qs, lq, alo);
            tmp[lq + 1] = mpn_addmul_1(tmp + 1, Qs, lq, ahi);
            lt = lq + 2;
            while (lt > 0 && tmp[lt - 1] == 0)
                lt--;
        }
        if (ls > 0)
        {
            cy = mpn_sub(tmp, tmp, lt, S, ls);
            FLINT_ASSERT(cy == 0);
            while (lt > 0 && tmp[lt - 1] == 0)
                lt--;
        }
        FLINT_ASSERT(lt > 0);

        if (k == 0)
        {
            /* p(0) = 1: S(0) = tmp (positive); q(0) = q0f[0] q0f[1]
               scales Qs -- and thereby Q(0,N) -- without touching
               T (the k = 0 term of the sum compensates) */
            flint_mpn_copyi(S, tmp, lt);
            ls = lt;
            if (q0f[0] != 1)
            {
                cy = mpn_mul_1(Qs, Qs, lq, q0f[0]);
                Qs[lq] = cy;
                lq += (cy != 0);
            }
            if (q0f[1] != 1)
            {
                cy = mpn_mul_1(Qs, Qs, lq, q0f[1]);
                Qs[lq] = cy;
                lq += (cy != 0);
            }
            FLINT_ASSERT(lq <= PI_BS_ALLOC);
            break;
        }
        else
        {
            ulong g1 = 6 * (ulong) k - 5;
            ulong g2 = 2 * (ulong) k - 1;
            ulong g3 = 6 * (ulong) k - 1;

#define PI_MUL1(D_, LD_, F_)                                    \
            do {                                                \
                cy = mpn_mul_1(D_, D_, LD_, F_);                \
                (D_)[LD_] = cy;                                 \
                LD_ += (cy != 0);                               \
            } while (0)

            /* S = g(k) tmp (sign: negative, tracked implicitly) */
#if FLINT_BITS == 64
            if (k <= 503000)   /* g = g1 g2 g3 < 2^64 */
            {
                ulong g = g1 * g2 * g3;
                cy = mpn_mul_1(S, tmp, lt, g);
                S[lt] = cy;
                ls = lt + (cy != 0);
            }
            else
#endif
            {
                cy = mpn_mul_1(S, tmp, lt, g1);
                S[lt] = cy;
                ls = lt + (cy != 0);
                PI_MUL1(S, ls, g2);
                PI_MUL1(S, ls, g3);
            }

            /* Qs *= k^3 C, one whole-limb factor at a time */
#if FLINT_BITS == 64
            if ((ulong) k <= UWORD(2642245))   /* k^3 < 2^64 */
            {
                ulong k3 = (ulong) k * (ulong) k * (ulong) k;
                PI_MUL1(Qs, lq, k3);
                PI_MUL1(Qs, lq, PI_C1 * PI_C2);
            }
            else
#endif
            {
                PI_MUL1(Qs, lq, (ulong) k);
                PI_MUL1(Qs, lq, (ulong) k);
                PI_MUL1(Qs, lq, (ulong) k);
                PI_MUL1(Qs, lq, PI_C1);
                PI_MUL1(Qs, lq, PI_C2);
            }

            /* G *= g(k) */
            if (need_p)
            {
#if FLINT_BITS == 64
                if (k <= 503000)
                {
                    ulong g = g1 * g2 * g3;
                    PI_MUL1(G, lg, g);
                }
                else
#endif
                {
                    PI_MUL1(G, lg, g1);
                    PI_MUL1(G, lg, g2);
                    PI_MUL1(G, lg, g3);
                }
            }
#undef PI_MUL1

            FLINT_ASSERT(lq + 3 <= PI_BS_ALLOC && ls <= PI_BS_ALLOC
                && lg <= PI_BS_ALLOC);
        }
    }

    /* P = (-1)^(number of k >= 1 factors) G */
    if (need_p)
        _fball_set_mpn(P, G, lg, (int) ((b - FLINT_MAX(a, 1)) & 1));
    _fball_set_mpn(Q, Qs, lq, 0);
    /* T = S(a): positive iff the block starts at k = 0 */
    _fball_set_mpn(T, S, ls, a != 0);
}

/* P, Q, T over [a, b) at precision n.  P is only computed when
   need_p is set: the parent's T merge uses the P of its LEFT child,
   while the P of a right child only feeds the parent's own P product,
   so the entire right spine of the tree skips P. */
static void
pi_bsplit(fball_t P, fball_t Q, fball_t T, slong a, slong b, slong n,
    int need_p, const ulong * q0f)
{
    if (b - a <= PI_BS_BLK)
    {
        pi_bsplit_basecase(P, Q, T, a, b, need_p, q0f);
    }
    else
    {
        slong m = a + (b - a) / 2;
        fball_t P2, Q2, T2;

        fball_init(P2);
        fball_init(Q2);
        fball_init(T2);

        pi_bsplit(P, Q, T, a, m, n, 1, q0f);
        pi_bsplit(P2, Q2, T2, m, b, n, need_p, q0f);

        fball_mul(T, T, Q2, n);         /* T1 Q2 */
        fball_mul(T2, T2, P, n);        /* P1 T2 */
        fball_add(T, T, T2, n);

        fball_mul(Q, Q, Q2, n);
        if (need_p)
            fball_mul(P, P, P2, n);

        fball_clear(P2);
        fball_clear(Q2);
        fball_clear(T2);
    }
}

/* shared driver: pi (or pi/4, with the extra 1/4 folded into the
   baked q(0) factor) = Q rsqrt(640320) / T */
static void
_fball_const_pi(fball_t pi, slong n, const ulong * q0f)
{
    /* terms: the tail after N terms is below
       82 (N+1) 2^(-47 N) relative to S; N = (64 n + 64)/47 + 3
       leaves a relative tail below 2^-(64 n + 100) */
    slong N = (FLINT_BITS * n + 64) / PI_BITS_PER_TERM + 3;
    slong nb = n;
    fball_t P, Q, T;

    fball_init(P);
    fball_init(Q);
    fball_init(T);

    pi_bsplit(P, Q, T, 0, N, nb, 0, q0f);

    /* the D (resp. D/4) normalization rode in on q(0) */
    fball_rsqrt_ui(P, 640320, n);
    fball_mul(Q, Q, P, n);
    fball_div(pi, Q, T, n);

    /* series tail */
    fball_add_error_2exp_rel(pi, -(FLINT_BITS * n + 96));

    fball_clear(P);
    fball_clear(Q);
    fball_clear(T);
}

void
fball_const_pi_chudnovsky(fball_t pi, slong n)
{
#if FLINT_BITS == 64
    const ulong q0f[2] = { PI_D1 * PI_D2, UWORD(1) };
#else
    const ulong q0f[2] = { PI_D1, PI_D2 };
#endif
    _fball_const_pi(pi, n, q0f);
}

/* pi/4 directly: D/4 = PI_D4_1 PI_D2 baked into q(0) */
static void
_fball_const_pi_div_4(fball_t pi, slong n)
{
#if FLINT_BITS == 64
    const ulong q0f[2] = { PI_D4_1 * PI_D2, UWORD(1) };
#else
    const ulong q0f[2] = { PI_D4_1, PI_D2 };
#endif
    _fball_const_pi(pi, n, q0f);
}

/* fixed_const_pi_div_4: verified floor truncations of pi/4 from a
   per-thread cache, extended on demand.  The cache stores
   floor(pi/4 B^nc) for the largest nc requested so far (rounded up
   geometrically); floors nest, so any shorter request is the top
   limbs of the cached entry.  Each extension recomputes through the
   Chudnovsky splitting with increasing guard limbs until the ball's
   radius determines the floor uniquely (_fixed_tab_store_floor). */

static FLINT_TLS_PREFIX nn_ptr _fixed_pi4_cache = NULL;
static FLINT_TLS_PREFIX slong _fixed_pi4_n = 0;
static FLINT_TLS_PREFIX int _fixed_pi4_cleanup_registered = 0;

void
_fixed_const_pi_div_4_clear(void)
{
    flint_free(_fixed_pi4_cache);
    _fixed_pi4_cache = NULL;
    _fixed_pi4_n = 0;
}

static void
_fixed_pi4_cleanup(void)
{
    _fixed_const_pi_div_4_clear();
    _fixed_pi4_cleanup_registered = 0;
}

void
fixed_const_pi_div_4(nn_ptr y, slong n)
{
    FLINT_ASSERT(n >= 1);

    if (n > _fixed_pi4_n)
    {
        slong nc = FLINT_MAX(FLINT_MAX(n, 2 * _fixed_pi4_n), 16);
        slong guard = 3;
        nn_ptr e = flint_malloc(nc * sizeof(ulong));
        int ok;

        for (;;)
        {
            fball_t v;

            fball_init(v);
            _fball_const_pi_div_4(v, nc + guard);
            ok = fball_get_fixed_floor(e, nc, v);
            fball_clear(v);
            if (ok)
                break;
            guard += 2 + guard / 2;
        }

        flint_free(_fixed_pi4_cache);
        _fixed_pi4_cache = e;
        _fixed_pi4_n = nc;

        if (!_fixed_pi4_cleanup_registered)
        {
            flint_register_cleanup_function(_fixed_pi4_cleanup);
            _fixed_pi4_cleanup_registered = 1;
        }
    }

    flint_mpn_copyi(y, _fixed_pi4_cache + (_fixed_pi4_n - n), n);
}
