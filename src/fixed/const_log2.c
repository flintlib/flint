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

/* Same hypergeometric series as arb_const_log2 -- this is Zuniga's
   d = 2 Ramanujan-type identity [Zun2025, Eq. 18], the fastest
   known single series for log 2 (binary splitting cost
   C_s = -4d/log|rho| = 0.9679 at rho = 1/3888):

       log(2) = (1/2160) sum_{k>=0} A(k) prod_{j=1}^{k} p(j)/q(j),

       A(k) = 1497 + 1794 k,
       p(k) = k (2k - 1),
       q(k) = 7776 k^2 + 7776 k + 1080 = 216 (6k+1)(6k+5),
       p(0) = q(0) = 1.

   The d = 4 and d = 6 series of the same paper ([Zun2025],
   Eqs. 21, 22, 23, 25; costs 1.13-1.30) were checked in this
   splitting paradigm and measured 15-42% slower across
   10^2..10^5 limbs, converging to the asymptotic cost ratio
   1.17, so none of them displaces Eq. 18.

   All terms positive; p(j)/q(j) < 1/3888 = 2^-11.9248 strictly, so
   the tail after N terms is below 2 (1 + 1.2 N) 2^(-11.9248 N)
   relatively.  q(k) fits a limb for k <= 4.7e7 (about 5.6e8 bits),
   and Q(a,b) grows at 12.9 + 2 log2 k bits per term against the
   11.92 bits/term decay, so the exact tree tops out at ~4x the
   target precision at 1e7 bits - the precision cap works harder
   here than for pi.

   The bsplit structure (P, Q, T; blocks of exact mpn integers with
   backward recurrence; need_p right-spine skip; the reconstruction
   scalar baked into q(0) at the leftmost leaf, which scales Q
   without touching T) is identical to const_pi.c, only the
   basecase constants and the positive sign pattern differ: S(k) = p(k) (A(k) Qs + S(k+1)) accumulates
   with an mpn_add instead of a subtraction. */

#ifndef L2_BS_BLK
#define L2_BS_BLK 64
#endif
#define L2_BS_ALLOC ((L2_BS_BLK * 80) / FLINT_BITS + 12)

#define L2_A0 UWORD(1497)
#define L2_A1 UWORD(1794)

/* 1/11.92 slightly overestimates 1/11.9248: safe term count */
#define L2_TERMS(n) (((FLINT_BITS * (n) + 96) * 100) / 1192 + 4)

static void
_fball_set_mpn2(fball_t x, nn_srcptr d, slong len, int negative)
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

/* exact P, Q, T over [a, b), 1 <= b - a <= L2_BS_BLK */
static void
log2_bsplit_basecase(fball_t P, fball_t Q, fball_t T, slong a, slong b,
    int need_p)
{
    ulong Qs[L2_BS_ALLOC], S[L2_BS_ALLOC], G[L2_BS_ALLOC],
          tmp[L2_BS_ALLOC];
    slong lq, ls, lg, lt, k;

    /* each linear factor of p(k) = k (2k-1) and
       q(k) = 216 (6k+1)(6k+5) must fit a limb; on 64-bit machines
       the fused products are used while they fit */
    FLINT_ASSERT((ulong) b <= (UWORD_MAX - 5) / 6);

    Qs[0] = 1; lq = 1;
    G[0] = 1; lg = 1;
    ls = 0;

    for (k = b - 1; k >= a; k--)
    {
        ulong cy;

        /* tmp = A(k) Qs + S; A(k) = A0 + A1 k needs two limbs once
           A1 k overflows */
        if ((ulong) k <= (UWORD_MAX - L2_A0) / L2_A1)
        {
            ulong A = L2_A0 + L2_A1 * (ulong) k;
            cy = mpn_mul_1(tmp, Qs, lq, A);
            tmp[lq] = cy;
            lt = lq + (cy != 0);
        }
        else
        {
            ulong ahi, alo;
            umul_ppmm(ahi, alo, L2_A1, (ulong) k);
            add_ssaaaa(ahi, alo, ahi, alo, UWORD(0), L2_A0);
            tmp[lq] = mpn_mul_1(tmp, Qs, lq, alo);
            tmp[lq + 1] = mpn_addmul_1(tmp + 1, Qs, lq, ahi);
            lt = lq + 2;
            while (lt > 0 && tmp[lt - 1] == 0)
                lt--;
        }
        if (ls > 0)
        {
            FLINT_ASSERT(ls <= lt);
            cy = mpn_add(tmp, tmp, lt, S, ls);
            tmp[lt] = cy;
            lt += (cy != 0);
        }

        if (k == 0)
        {
            /* p(0) = 1; q(0) = 2160 bakes the reconstruction
               denominator log(2) = T / (2160 Q) into Q itself */
            flint_mpn_copyi(S, tmp, lt);
            ls = lt;
            cy = mpn_mul_1(Qs, Qs, lq, UWORD(2160));
            Qs[lq] = cy;
            lq += (cy != 0);
            FLINT_ASSERT(lq <= L2_BS_ALLOC);
            break;
        }
        else
        {
#define L2_MUL1(D_, LD_, F_)                                    \
            do {                                                \
                cy = mpn_mul_1(D_, D_, LD_, F_);                \
                (D_)[LD_] = cy;                                 \
                LD_ += (cy != 0);                               \
            } while (0)

            /* S = p(k) tmp = k (2k-1) tmp */
#if FLINT_BITS == 64
            if ((ulong) k <= UWORD(3037000499))   /* p < 2^64 */
            {
                ulong p = (ulong) k * (2 * (ulong) k - 1);
                cy = mpn_mul_1(S, tmp, lt, p);
                S[lt] = cy;
                ls = lt + (cy != 0);
            }
            else
#endif
            {
                cy = mpn_mul_1(S, tmp, lt, (ulong) k);
                S[lt] = cy;
                ls = lt + (cy != 0);
                L2_MUL1(S, ls, 2 * (ulong) k - 1);
            }

            /* Qs *= q(k) = 216 (6k+1)(6k+5) */
#if FLINT_BITS == 64
            if ((ulong) k <= UWORD(46000000))   /* q < 2^64 */
            {
                ulong q = 216 * ((6 * (ulong) k + 1)
                    * (6 * (ulong) k + 5));
                L2_MUL1(Qs, lq, q);
            }
            else
#endif
            {
                L2_MUL1(Qs, lq, 6 * (ulong) k + 1);
                L2_MUL1(Qs, lq, 6 * (ulong) k + 5);
                L2_MUL1(Qs, lq, UWORD(216));
            }

            /* G *= p(k) */
            if (need_p)
            {
#if FLINT_BITS == 64
                if ((ulong) k <= UWORD(3037000499))
                {
                    ulong p = (ulong) k * (2 * (ulong) k - 1);
                    L2_MUL1(G, lg, p);
                }
                else
#endif
                {
                    L2_MUL1(G, lg, (ulong) k);
                    L2_MUL1(G, lg, 2 * (ulong) k - 1);
                }
            }
#undef L2_MUL1

            FLINT_ASSERT(lq + 3 <= L2_BS_ALLOC && ls <= L2_BS_ALLOC
                && lg <= L2_BS_ALLOC);
        }
    }

    if (need_p)
        _fball_set_mpn2(P, G, lg, 0);
    _fball_set_mpn2(Q, Qs, lq, 0);
    _fball_set_mpn2(T, S, ls, 0);
}

static void
log2_bsplit(fball_t P, fball_t Q, fball_t T, slong a, slong b, slong n,
    int need_p)
{
    if (b - a <= L2_BS_BLK)
    {
        log2_bsplit_basecase(P, Q, T, a, b, need_p);
    }
    else
    {
        slong m = a + (b - a) / 2;
        fball_t P2, Q2, T2;

        fball_init(P2);
        fball_init(Q2);
        fball_init(T2);

        log2_bsplit(P, Q, T, a, m, n, 1);
        log2_bsplit(P2, Q2, T2, m, b, n, need_p);

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

void
fball_const_log2(fball_t res, slong n)
{
    slong N = L2_TERMS(n);
    fball_t P, Q, T;

    fball_init(P);
    fball_init(Q);
    fball_init(T);

    log2_bsplit(P, Q, T, 0, N, n, 0);

    /* the 2160 reconstruction denominator rode in on q(0) */
    fball_div(res, T, Q, n);

    /* series tail */
    fball_add_error_2exp_rel(res, -(FLINT_BITS * n + 96));

    fball_clear(P);
    fball_clear(Q);
    fball_clear(T);
}

/* fixed_const_log2: verified floor truncations of log(2) from a
   per-thread cache, exactly as fixed_const_pi_div_4 (see
   const_pi.c). */

static FLINT_TLS_PREFIX nn_ptr _fixed_log2_cache = NULL;
static FLINT_TLS_PREFIX slong _fixed_log2_n = 0;
static FLINT_TLS_PREFIX int _fixed_log2_cleanup_registered = 0;

void
_fixed_const_log2_clear(void)
{
    flint_free(_fixed_log2_cache);
    _fixed_log2_cache = NULL;
    _fixed_log2_n = 0;
}

static void
_fixed_log2_cleanup(void)
{
    _fixed_const_log2_clear();
    _fixed_log2_cleanup_registered = 0;
}

void
fixed_const_log2(nn_ptr y, slong n)
{
    FLINT_ASSERT(n >= 1);

    if (n > _fixed_log2_n)
    {
        slong nc = FLINT_MAX(FLINT_MAX(n, 2 * _fixed_log2_n), 16);
        slong guard = 3;
        nn_ptr e = flint_malloc(nc * sizeof(ulong));
        int ok;

        for (;;)
        {
            fball_t v;

            fball_init(v);
            fball_const_log2(v, nc + guard);
            ok = fball_get_fixed_floor(e, nc, v);
            fball_clear(v);
            if (ok)
                break;
            guard += 2 + guard / 2;
        }

        flint_free(_fixed_log2_cache);
        _fixed_log2_cache = e;
        _fixed_log2_n = nc;

        if (!_fixed_log2_cleanup_registered)
        {
            flint_register_cleanup_function(_fixed_log2_cleanup);
            _fixed_log2_cleanup_registered = 1;
        }
    }

    flint_mpn_copyi(y, _fixed_log2_cache + (_fixed_log2_n - n), n);
}
