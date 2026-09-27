/*
    Copyright (C) 2012, 2013, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "mp_real.h"
#include "impl.h"
#include "arb_hypgeom.h"

/* Euler's constant by the Brent-McMillan algorithm (the refinement B3
   with the asymptotic series for I_0 K_0), as in arb_const_euler, in
   mp_real arithmetic:

       gamma = S_0 / I_0 - K_0 / I_0^2 - log m + O(e^{-8m}),

       I_0 = sum_{k>=0} (m^k / k!)^2,   S_0 = sum_{k>=0} (m^k / k!)^2 H_k,
       K_0 = I_0(2m) K_0(2m) = (1/(4m)) sum_{k<2m} (-1)^k ((2k)!)^3
                               / ((k!)^4 (16 m)^(2k)),

   the error bounded by 24 e^{-8m} [BJ2013].  The first two sums are
   evaluated by one binary splitting over N = ceil(4.970626 m) + 1 terms
   in dual numbers (see below; arb's euler_bsplit_1 uses the quantities
   P, Q, T, C, D, V instead, with about twice the multiplication
   cost), the third by a plain P, Q, T splitting over 2m terms at half
   the precision.  Exact fmpz products
   below EULER_LEAF terms, mp_real_t balls (exact until they outgrow the working
   precision) above.

   THE LOGARITHM.  m only needs m >= 0.0866434 b + 1 for b bits
   (8m >= b log 2); it is rounded up to the next m that is smooth over
   a small set of primes so that log m is a combination of a few
   logarithms, each from a fast series.  The sets tried (log m from
   Zuniga's series [Zun2025] of the listed ratios, combined exactly):

       {2}         log 2 (mp_real_const_log2)
       {2,3}       9/8, 256/243
       {2,3,5}     81/80, 32805/32768, 25/24
       {2,3,5,7}   2401/2400, 4375/4374, 225/224, 64/63

   The logarithms cost 6-7% of the total for any set; what matters more
   is m itself, through the size of the odd part of m (see
   _euler_cost), so every m in [m0, 1.5 m0] smooth over {2,3,5,7} is
   scored and the cheapest taken (over 10^4 - 3 10^6 bits, 0.2% above
   the best fixed set at each precision on geometric average, against
   2-3% for any one of {2,3}, {2,3,5}, {2,3,5,7} and 26% for {2}). */

#define EULER_LEAF 16

/* right children above this many limbs are freed after the merge */
#define EULER_KEEP 2048

/* ---- mp_real helpers ---- */

static void
_mp_real_set_fmpz(mp_real_t x, const fmpz_t f)
{
    slong n = fmpz_size(f);
    fmpz_t a;
    nn_ptr d;

    if (n == 0)
    {
        mp_real_zero(x);
        return;
    }
    if (!COEFF_IS_MPZ(*f))
    {
        ulong u = FLINT_UABS(*f);
        _mp_real_set_mpn_2exp(x, &u, 1, 0);
        if (*f < 0)
            mp_real_neg(x, x);
        return;
    }
    fmpz_init(a);
    fmpz_abs(a, f);
    d = flint_malloc(n * sizeof(ulong));
    fmpz_get_ui_array(d, n, a);
    _mp_real_set_mpn_2exp(x, d, n, 0);
    if (fmpz_sgn(f) < 0)
        mp_real_neg(x, x);
    flint_free(d);
    fmpz_clear(a);
}

static mp_real_struct *
_tmp_init(slong k)
{
    slong i;
    mp_real_struct * tmp = flint_malloc(k * sizeof(mp_real_struct));
    for (i = 0; i < k; i++)
        mp_real_init(tmp + i);
    return tmp;
}

static void
_tmp_clear(mp_real_struct * tmp, slong k)
{
    slong i;
    for (i = 0; i < k; i++)
        mp_real_clear(tmp + i);
    flint_free(tmp);
}

/* ---- S_0 and I_0 ---- */

/* With F(e) = sum_{k>=0} m^(2k) / prod_{j=1}^k (j + e)^2 over the dual
   numbers (e^2 = 0), I_0 = F(0) and S_0 = -F'(0) / 2, so a single P, Q, T
   splitting of F over Z[e] / (e^2) gives both sums.  Q(e) = D(e)^2 with
   D(e) = prod (k + 1 + e) of half the size, so the state is P (a
   scalar), D = D0 + D1 e and T = T0 + T1 e; the merge

       T = P1 T2 + Q2 T1,  D = D1 D2,  P = P1 P2,  Q2 = D2^2

   costs 5 full-size and 5 half-size products (against 12 mostly
   full-size ones for arb's P, Q, T, C, D, V).  At the end, with
   X = T0 + D0^2,

       I_0 = X / D0^2,   S_0 / I_0 = (2 T0 D1 - T1 D0) / (2 D0 X).

   T1 >= 0 (the e-part of each partial sum T / Q is negative, but Q1 / Q0
   dominates it), so no merge loses more than a few bits to
   cancellation.  The final divisions cost about 1.5% of the total;
   fusing them into one saves a half-precision division but needs more
   full-precision products, so they are kept apart. */

typedef struct
{
    mp_real_struct P, D0, D1, T0, T1;
}
e1_struct;

typedef struct
{
    fmpz_t P, D0, D1, T0, T1;
}
e1z_struct;

static void
e1z_init(e1z_struct * s)
{
    fmpz_init(s->P); fmpz_init(s->D0); fmpz_init(s->D1);
    fmpz_init(s->T0); fmpz_init(s->T1);
}

static void
e1z_clear(e1z_struct * s)
{
    fmpz_clear(s->P); fmpz_clear(s->D0); fmpz_clear(s->D1);
    fmpz_clear(s->T0); fmpz_clear(s->T1);
}

static void
e1_init(e1_struct * s)
{
    mp_real_init(&s->P); mp_real_init(&s->D0); mp_real_init(&s->D1);
    mp_real_init(&s->T0); mp_real_init(&s->T1);
}

static void
e1_clear(e1_struct * s)
{
    mp_real_clear(&s->P); mp_real_clear(&s->D0); mp_real_clear(&s->D1);
    mp_real_clear(&s->T0); mp_real_clear(&s->T1);
}

/* exact values over [a, b): single terms P = T = m^2, D = k + 1 + e */
static void
e1z_bsplit(e1z_struct * s, slong a, slong b, ulong m)
{
    if (b - a == 1)
    {
        fmpz_set_ui(s->P, m);
        fmpz_mul_ui(s->P, s->P, m);
        fmpz_set(s->T0, s->P);
        fmpz_zero(s->T1);
        fmpz_set_ui(s->D0, a + 1);
        fmpz_one(s->D1);
    }
    else
    {
        e1z_struct R;
        fmpz_t q0, q1, t;
        slong mid = a + (b - a) / 2;

        e1z_init(&R);
        fmpz_init(q0);
        fmpz_init(q1);
        fmpz_init(t);
        e1z_bsplit(s, a, mid, m);
        e1z_bsplit(&R, mid, b, m);

        /* Q2 = D2^2 */
        fmpz_mul(q0, R.D0, R.D0);
        fmpz_mul(q1, R.D0, R.D1);
        fmpz_mul_2exp(q1, q1, 1);

        /* T1 = P1 T2_1 + Q2_0 T1_1 + Q2_1 T1_0 */
        fmpz_mul(s->T1, s->T1, q0);
        fmpz_addmul(s->T1, q1, s->T0);
        fmpz_addmul(s->T1, s->P, R.T1);
        /* T0 = P1 T2_0 + Q2_0 T1_0 */
        fmpz_mul(s->T0, s->T0, q0);
        fmpz_addmul(s->T0, s->P, R.T0);

        /* D = D1 D2 */
        fmpz_mul(t, s->D0, R.D1);
        fmpz_mul(s->D1, s->D1, R.D0);
        fmpz_add(s->D1, s->D1, t);
        fmpz_mul(s->D0, s->D0, R.D0);

        fmpz_mul(s->P, s->P, R.P);

        e1z_clear(&R);
        fmpz_clear(q0);
        fmpz_clear(q1);
        fmpz_clear(t);
    }
}

/* subtrees of at least this many terms run their halves on two
   threads when the thread budget allows; the right one then gets its
   own level temporaries and scratch, freed when it is done */
#define EULER_PAR_TERMS 4096

static void e1_bsplit(e1_struct * s, slong a, slong b, ulong m, int need,
    slong wp, mp_real_struct * tmp, mp_real_struct * sc);
static void e2_bsplit(mp_real_t P, mp_real_t Q, mp_real_t T, slong a, slong b,
    ulong m, int need, slong wp, mp_real_struct * tmp, mp_real_struct * u);

typedef struct
{
    mp_real_struct * s;       /* e1: an e1_struct; e2: P, Q, T */
    mp_real_struct * tmp, * sc;
    slong a, b, wp;
    ulong m;
    int need, e2;
}
euler_job;

static void
_euler_job(void * arg)
{
    euler_job * J = (euler_job *) arg;
    mp_real_struct * tmp = J->tmp, * sc = J->sc;
    slong k = 0, per = J->e2 ? 3 : 5, nsc = J->e2 ? 1 : 6;

    if (tmp == NULL)
    {
        k = nsc + per * (FLINT_BIT_COUNT(J->b - J->a) + 2);
        sc = _tmp_init(k);
        tmp = sc + nsc;
    }

    if (J->e2)
        e2_bsplit(J->s, J->s + 1, J->s + 2, J->a, J->b, J->m, J->need,
            J->wp, tmp, sc);
    else
        e1_bsplit((e1_struct *) J->s, J->a, J->b, J->m, J->need, J->wp,
            tmp, sc);

    if (k != 0)
        _tmp_clear(sc, k);
}

/* the rule of MP_REAL_PAR_CAP for [a, b): 0 serial, 1 fork the halves,
   2 above the cap (halves serial, parallel merge); the bits of T over
   [a, b) estimated as sum log2(m^2 (k+1)^2) for sum 1, sum log2(32 k^3
   m^2) for sum 2 */
static int
_euler_par(slong a, slong b, ulong m, slong wp, int e2)
{
    double est;

    if (b - a < EULER_PAR_TERMS || flint_get_num_threads() < 2)
        return 0;
    if (e2)
        est = (double) (b - a) * (3.0 * log2((double) b)
            + 2.0 * log2((double) m) + 5.0);
    else
        est = (double) (b - a) * 2.0 * log2((double) m * (double) b);
    return (est > MP_REAL_PAR_CAP * FLINT_BITS * (double) wp) ? 2 : 1;
}

/* the halves [a, mid), [mid, b) into s and R (1, need) on two threads */
static void
_euler_fork(mp_real_struct * s, mp_real_struct * R, slong a, slong mid,
    slong b, ulong m, int need, slong wp, mp_real_struct * tmp,
    mp_real_struct * sc, int e2)
{
    euler_job L, Rj;

    L.s = s; L.tmp = tmp; L.sc = sc; L.a = a; L.b = mid; L.wp = wp;
    L.m = m; L.need = 1; L.e2 = e2;
    Rj.s = R; Rj.tmp = NULL; Rj.sc = NULL; Rj.a = mid; Rj.b = b;
    Rj.wp = wp; Rj.m = m; Rj.need = need; Rj.e2 = e2;
    _mp_real_parallel_pair(_euler_job, &L, _euler_job, &Rj);
}

/* the sum 1 merge (see e1_bsplit) on two threads, in two rounds:
   {Q2_0 = D2_0^2, u = P1 T2_1}, {Q2_1 = 2 D2_0 D2_1, w = P1 T2_0}, then
   {T1 = u + Q2_0 T1_1 + Q2_1 T1_0}, {w += Q2_0 T1_0, D = D1 D2}; sc
   holds 6 scratch variables */
typedef struct
{
    e1_struct * s, * R;
    mp_real_struct * sc;
    slong sh, wp;
    int round;
}
e1_merge_struct;

static void
_e1_merge_x(void * arg)
{
    e1_merge_struct * M = (e1_merge_struct *) arg;
    e1_struct * s = M->s, * R = M->R;
    mp_real_struct * q0 = M->sc, * q1 = M->sc + 1, * u = M->sc + 2,
        * v = M->sc + 3;
    slong wp = M->wp;

    if (M->round == 0)
    {
        mp_real_mul(q0, &R->D0, &R->D0, wp);
        mp_real_mul(u, &s->P, &R->T1, wp);
        mp_real_mul_2exp_si(u, u, M->sh);
    }
    else
    {
        mp_real_mul(v, q0, &s->T1, wp);
        mp_real_add(u, u, v, wp);
        mp_real_mul(v, q1, &s->T0, wp);
        mp_real_add(&s->T1, u, v, wp);
    }
}

static void
_e1_merge_y(void * arg)
{
    e1_merge_struct * M = (e1_merge_struct *) arg;
    e1_struct * s = M->s, * R = M->R;
    mp_real_struct * q0 = M->sc, * q1 = M->sc + 1, * w = M->sc + 4,
        * v = M->sc + 5;
    slong wp = M->wp;

    if (M->round == 0)
    {
        mp_real_mul(q1, &R->D0, &R->D1, wp);
        mp_real_mul_2exp_si(q1, q1, 1);
        mp_real_mul(w, &s->P, &R->T0, wp);
        mp_real_mul_2exp_si(w, w, M->sh);
    }
    else
    {
        mp_real_mul(v, q0, &s->T0, wp);
        mp_real_add(w, w, v, wp);
        mp_real_mul(v, &s->D0, &R->D1, wp);
        mp_real_mul(&s->D1, &s->D1, &R->D0, wp);
        mp_real_add(&s->D1, &s->D1, v, wp);
        mp_real_mul(&s->D0, &s->D0, &R->D0, wp);
    }
}

static void
_e1_merge_par(e1_struct * s, e1_struct * R, int need, slong sh, slong wp,
    mp_real_struct * sc)
{
    e1_merge_struct M;

    M.s = s; M.R = R; M.sc = sc; M.sh = sh; M.wp = wp;
    M.round = 0;
    _mp_real_parallel_pair(_e1_merge_x, &M, _e1_merge_y, &M);
    M.round = 1;
    _mp_real_parallel_pair(_e1_merge_x, &M, _e1_merge_y, &M);
    mp_real_swap(&s->T0, sc + 4);
    if (need)
        mp_real_mul(&s->P, &s->P, &R->P, wp);
}

/* the same in mp_real_t balls at wp limbs; P only when need (not on the right
   spine).  tmp: the right child, 5 per level below; sc: 6 scratch
   variables shared by all levels (one merge runs at a time) */
static void
e1_bsplit(e1_struct * s, slong a, slong b, ulong m, int need, slong wp,
    mp_real_struct * tmp, mp_real_struct * sc)
{
    if (b - a <= EULER_LEAF)
    {
        e1z_struct z;
        e1z_init(&z);
        e1z_bsplit(&z, a, b, m);
        /* P = m^(2 (b - a)) is stored without its power of two */
        fmpz_fdiv_q_2exp(z.P, z.P, 2 * flint_ctz(m) * (b - a));
        _mp_real_set_fmpz(&s->P, z.P);
        _mp_real_set_fmpz(&s->D0, z.D0);
        _mp_real_set_fmpz(&s->D1, z.D1);
        _mp_real_set_fmpz(&s->T0, z.T0);
        _mp_real_set_fmpz(&s->T1, z.T1);
        e1z_clear(&z);
    }
    else
    {
        e1_struct * R = (e1_struct *) tmp;
        mp_real_struct * q0 = sc, * q1 = sc + 1, * u = sc + 2, * v = sc + 3;
        slong mid = a + (b - a) / 2;

        slong sh = 2 * flint_ctz(m) * (mid - a);   /* P1's power of two */
        int par = _euler_par(a, b, m, wp, 0);

        if (par == 1)
            _euler_fork((mp_real_struct *) s, (mp_real_struct *) R, a, mid, b,
                m, need, wp, tmp + 5, sc, 0);
        else
        {
            e1_bsplit(s, a, mid, m, 1, wp, tmp + 5, sc);
            e1_bsplit(R, mid, b, m, need, wp, tmp + 5, sc);
        }

        if (par == 2)
        {
            _e1_merge_par(s, R, need, sh, wp, sc);
            goto release;
        }

        /* Q2 = D2^2 */
        mp_real_mul(q0, &R->D0, &R->D0, wp);
        mp_real_mul(q1, &R->D0, &R->D1, wp);
        mp_real_mul_2exp_si(q1, q1, 1);

        /* T1 = P1 T2_1 + Q2_0 T1_1 + Q2_1 T1_0 */
        mp_real_mul(u, &s->P, &R->T1, wp);
        mp_real_mul_2exp_si(u, u, sh);
        mp_real_mul(v, q0, &s->T1, wp);
        mp_real_add(u, u, v, wp);
        mp_real_mul(v, q1, &s->T0, wp);
        mp_real_add(&s->T1, u, v, wp);

        /* T0 = P1 T2_0 + Q2_0 T1_0 */
        mp_real_mul(u, &s->P, &R->T0, wp);
        mp_real_mul_2exp_si(u, u, sh);
        mp_real_mul(v, q0, &s->T0, wp);
        mp_real_add(&s->T0, u, v, wp);

        /* D = D1 D2 */
        mp_real_mul(u, &s->D0, &R->D1, wp);
        mp_real_mul(v, &s->D1, &R->D0, wp);
        mp_real_add(&s->D1, u, v, wp);
        /* in place: swapping with the shared scratch would move large
           buffers down into the deep levels, where they stay */
        mp_real_mul(&s->D0, &s->D0, &R->D0, wp);

        if (need)
            mp_real_mul(&s->P, &s->P, &R->P, wp);

release:
        /* free a large right child: this level's buffers are otherwise
           held through all the merges above */
        if (R->T0.alloc > EULER_KEEP)
        {
            e1_clear(R);
            e1_init(R);
        }
    }
}

/* ---- K_0 ---- */

/* exact P, Q, T over [a, b): p(0) = 1, q(0) = 4m, p(k) = (2k-1)^3,
   q(k) = 32 k m^2 */
static void
e2z_bsplit(fmpz_t P, fmpz_t Q, fmpz_t T, slong a, slong b, ulong m)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            fmpz_one(P);
            fmpz_set_ui(Q, 4 * m);
        }
        else
        {
            fmpz_set_ui(P, 2 * a - 1);
            fmpz_pow_ui(P, P, 3);
            fmpz_set_ui(Q, 32 * (ulong) a);
            fmpz_mul_ui(Q, Q, m);
            fmpz_mul_ui(Q, Q, m);
        }
        fmpz_set(T, P);
    }
    else
    {
        fmpz_t P2, Q2, T2;
        slong mid = a + (b - a) / 2;

        fmpz_init(P2);
        fmpz_init(Q2);
        fmpz_init(T2);
        e2z_bsplit(P, Q, T, a, mid, m);
        e2z_bsplit(P2, Q2, T2, mid, b, m);
        fmpz_mul(T, T, Q2);
        fmpz_mul(T2, T2, P);
        fmpz_add(T, T, T2);
        fmpz_mul(P, P, P2);
        fmpz_mul(Q, Q, Q2);
        fmpz_clear(P2);
        fmpz_clear(Q2);
        fmpz_clear(T2);
    }
}

/* P, Q, T at wp limbs; P, Q, T of the left child must be consecutive
   (the caller's P + 1 == Q, P + 2 == T) for the threaded split */
static void
e2_bsplit(mp_real_t P, mp_real_t Q, mp_real_t T, slong a, slong b, ulong m,
    int need, slong wp, mp_real_struct * tmp, mp_real_struct * u)
{
    if (b - a <= EULER_LEAF)
    {
        fmpz_t zP, zQ, zT;
        fmpz_init(zP);
        fmpz_init(zQ);
        fmpz_init(zT);
        e2z_bsplit(zP, zQ, zT, a, b, m);
        _mp_real_set_fmpz(P, zP);
        _mp_real_set_fmpz(Q, zQ);
        _mp_real_set_fmpz(T, zT);
        fmpz_clear(zP);
        fmpz_clear(zQ);
        fmpz_clear(zT);
    }
    else
    {
        mp_real_struct * P2 = tmp, * Q2 = tmp + 1, * T2 = tmp + 2;
        slong mid = a + (b - a) / 2;

        int par = _euler_par(a, b, m, wp, 1);

        FLINT_ASSERT(Q == P + 1 && T == P + 2);
        if (par == 1)
            _euler_fork(P, P2, a, mid, b, m, need, wp, tmp + 3, u, 1);
        else
        {
            e2_bsplit(P, Q, T, a, mid, m, 1, wp, tmp + 3, u);
            e2_bsplit(P2, Q2, T2, mid, b, m, need, wp, tmp + 3, u);
        }

        _mp_real_pqt_merge(P, Q, T, P2, Q2, T2, need, wp, par == 2);

        if (T2->alloc > EULER_KEEP)
        {
            mp_real_clear(P2); mp_real_init(P2);
            mp_real_clear(Q2); mp_real_init(Q2);
            mp_real_clear(T2); mp_real_init(T2);
        }
    }
}

/* ---- log m ---- */

/* the sets: primes, ratios u_i / v_i, and den X^{-1} (row j: log p_j
   = (1/den) sum_i A_ji log(u_i/v_i)) */
typedef struct
{
    int np;
    ulong u[4], v[4];
    slong A[16];
    int den;
}
euler_logset_t;

static const euler_logset_t euler_logsets[4] = {
    { 1, { 2 }, { 1 }, { 1 }, 1 },
    { 2, { 9, 256 }, { 8, 243 }, { 5, 2, 8, 3 }, 1 },
    { 3, { 81, 32805, 25 }, { 80, 32768, 24 },
      { 17, -7, 12, 27, -11, 19, 39, -16, 28 }, 1 },
    { 4, { 2401, 4375, 225, 64 }, { 2400, 4374, 224, 63 },
      { 34, -10, 54, 72, 54, -16, 86, 114, 79, -23, 126, 167,
        96, -28, 152, 202 }, 2 },
};

static const ulong euler_primes[4] = { 2, 3, 5, 7 };

/* the estimated cost of the splittings for the parameter m: m
   (EULER_ALPHA + log2 of the odd part of m), as the power of two in
   m^2 is free (P is stored without it, and zero limbs of the exact
   products are stripped) while the rest costs about 1/EULER_ALPHA of
   the splitting per bit (measured at 10^6 bits over the 7-smooth m
   within 35% of the minimum: within 2%) */
#define EULER_ALPHA 61.0

static double
_euler_cost(ulong m)
{
    ulong o = m >> flint_ctz(m);

    /* for m a power of two P is free altogether (measured: 5% below
       the fit at odd part 3) */
    if (o == 1)
        return (double) m * (EULER_ALPHA - 3.0);
    return (double) m * (EULER_ALPHA + log2((double) o));
}

/* relative cost of log m for the sets {2}, {2,3}, {2,3,5}, {2,3,5,7}
   against the rest (measured at 10^5 - 10^7 bits) */
static const double euler_log_cost[4] = { 0.063, 0.066, 0.074, 0.072 };

/* the m >= m0 smooth over (a prefix of) 2, 3, 5, 7 minimizing the
   estimated cost; np = 1 ... 4 restricts to the first np primes, np = 0
   lets the set follow m */
static ulong
_euler_choose_m(int * set, ulong m0, int np)
{
    ulong k, t, best = 0;
    double c, bc = 0.0;
    int j, need;

    for (k = m0; k <= m0 + m0 / 2 + 2; k++)
    {
        t = k;
        need = 1;
        for (j = 0; j < (np ? np : 4); j++)
            while (t % euler_primes[j] == 0)
            {
                t /= euler_primes[j];
                need = j + 1;
            }
        if (t != 1)
            continue;
        if (np)
            need = np;
        c = _euler_cost(k) * (1.0 + euler_log_cost[need - 1]);
        if (best == 0 || c < bc)
        {
            best = k;
            bc = c;
            *set = need - 1;
        }
    }

    if (best == 0)   /* np = 1: a power of two */
    {
        best = UWORD(1) << FLINT_BIT_COUNT(m0 - 1);
        *set = 0;
    }
    return best;
}

typedef struct
{
    mp_real_struct * v;
    const euler_logset_t * L;
    slong idx[4];
    slong wp;
}
euler_log_work;

static void
_euler_log_worker(slong j, void * arg)
{
    euler_log_work * W = (euler_log_work *) arg;
    slong i = W->idx[j];
    fmpz_t u, v;

    if (W->L->np == 1)
    {
        _mp_real_const_log2_compute(W->v + i, W->wp);
        return;
    }
    fmpz_init_set_ui(u, W->L->u[i]);
    fmpz_init_set_ui(v, W->L->v[i]);
    _mp_real_log_ratio_zuniga(W->v + i, u, v, W->wp);
    fmpz_clear(u);
    fmpz_clear(v);
}

/* log m for m smooth over the primes of the set: the needed series on
   up to np threads, the costliest (smallest u) first */
static void
_euler_log(mp_real_t res, ulong m, const euler_logset_t * L, slong wp)
{
    slong e[4], w[4], i, j, cnt;
    mp_real_struct v[4];
    mp_real_t t;
    euler_log_work W;

    for (j = 0; j < L->np; j++)
    {
        e[j] = 0;
        while (m % euler_primes[j] == 0)
        {
            m /= euler_primes[j];
            e[j]++;
        }
    }
    FLINT_ASSERT(m == 1);

    /* w_i = sum_j e_j A_ji */
    cnt = 0;
    for (i = 0; i < L->np; i++)
    {
        w[i] = 0;
        for (j = 0; j < L->np; j++)
            w[i] += e[j] * L->A[j * L->np + i];
        mp_real_init(v + i);
        if (w[i] != 0)
            W.idx[cnt++] = i;
    }

    /* by increasing u */
    for (i = 1; i < cnt; i++)
        for (j = i; j > 0 && L->u[W.idx[j]] < L->u[W.idx[j - 1]]; j--)
            FLINT_SWAP(slong, W.idx[j], W.idx[j - 1]);

    W.v = v;
    W.L = L;
    W.wp = wp;
    _mp_real_parallel_tasks(_euler_log_worker, &W, cnt);

    mp_real_init(t);
    mp_real_zero(res);
    for (i = 0; i < L->np; i++)
    {
        if (w[i] != 0)
        {
            if (w[i] < 0)
                mp_real_submul_ui(res, res, v + i, FLINT_UABS(w[i]), wp);
            else
                mp_real_addmul_ui(res, res, v + i, FLINT_UABS(w[i]), wp);
        }
        mp_real_clear(v + i);
    }

    if (L->den == 2)
        mp_real_mul_2exp_si(res, res, -1);

    mp_real_clear(t);
}

/* ---- the driver ---- */

/* x += [-2^e, 2^e] */
static void
_mp_real_add_error_2exp(mp_real_t x, slong e)
{
    slong q = e >> (FLINT_BITS == 64 ? 6 : 5);
    _mp_real_add_error_ulps_at(x, (double) (UWORD(1) << (e - q * FLINT_BITS)), q);
}



/* The main sum runs first, with nothing else live; it leaves res and
   the two half-precision values Q0^2, X^2, and K_0 and log m (each
   with its own scratch, freed at the end) follow. */
void
_mp_real_const_euler_tune(mp_real_t res, slong n, int set)
{
    const euler_logset_t * LS;
    slong bits = FLINT_BITS * n + 10, N, M, wp, wp2, depth, k;
    ulong m;
    e1_struct S;
    mp_real_struct * tmp;
    mp_real_t A, B, t, u;

    m = (ulong) (0.086643397569993163677 * bits) + 1;
    m = _euler_choose_m(&set, m, set < 0 ? 0 : set + 1);
    LS = euler_logsets + set;

    N = (slong) ceil(4.970626 * (double) m) + 1;
    M = 2 * m;

    wp = n + (2 * FLINT_BIT_COUNT(m)) / FLINT_BITS + 1;
    wp2 = n / 2 + (2 * FLINT_BIT_COUNT(m)) / FLINT_BITS + 2;

    mp_real_init(A);
    mp_real_init(B);
    mp_real_init(t);
    mp_real_init(u);

    /* S_0, I_0: 5 per level and 6 shared scratch */
    depth = FLINT_BIT_COUNT(N) + 2;
    k = 5 * depth + 6;
    tmp = _tmp_init(k);
    e1_init(&S);
    e1_bsplit(&S, 0, N, m, 0, wp, tmp + 6, tmp);
    _tmp_clear(tmp, k);

    /* S_0 / I_0 = (2 T0 D1 - T1 D0) / (2 D0 X), X = T0 + D0^2 */
    mp_real_mul(t, &S.T0, &S.D1, wp);
    mp_real_mul_2exp_si(t, t, 1);
    mp_real_mul(u, &S.T1, &S.D0, wp);
    mp_real_sub(t, t, u, wp);
    mp_real_mul(u, &S.D0, &S.D0, wp);                 /* Q0 = D0^2 */
    mp_real_mul(A, u, u, wp2);                        /* Q0^2 */
    mp_real_add(&S.T0, &S.T0, u, wp);                 /* X */
    mp_real_mul(B, &S.T0, &S.T0, wp2);                /* X^2 */
    mp_real_mul(u, &S.T0, &S.D0, wp);
    mp_real_mul_2exp_si(u, u, 1);
    e1_clear(&S);
    mp_real_div(res, t, u, wp);
    mp_real_clear(t);
    mp_real_clear(u);
    mp_real_init(t);
    mp_real_init(u);

    /* K_0 / I_0^2 = Q0^2 T2 / (X^2 Q2) at half precision; 3 per level
       and one shared scratch */
    {
        mp_real_struct PQT[3];    /* consecutive, for the threaded split */

        mp_real_init(PQT);
        mp_real_init(PQT + 1);
        mp_real_init(PQT + 2);
        depth = FLINT_BIT_COUNT(M) + 2;
        k = 3 * depth + 1;
        tmp = _tmp_init(k);
        e2_bsplit(PQT, PQT + 1, PQT + 2, 0, M, m, 0, wp2, tmp + 1, tmp);
        _tmp_clear(tmp, k);
        mp_real_clear(PQT);
        mp_real_mul(A, A, PQT + 2, wp2);
        mp_real_mul(B, B, PQT + 1, wp2);
        mp_real_clear(PQT + 1);
        mp_real_clear(PQT + 2);
        mp_real_div(t, A, B, wp2);
        mp_real_sub(res, res, t, wp);
    }

    /* log m */
    _euler_log(t, m, LS, wp);
    mp_real_sub(res, res, t, wp);

    /* 24 e^{-8m} <= 2^(4.6 - 11.5415 m) */
    _mp_real_add_error_2exp(res, (slong) ceil(4.6 - 11.5415 * (double) m));

    mp_real_clear(A);
    mp_real_clear(B);
    mp_real_clear(t);
    mp_real_clear(u);
}

FLINT_DLL extern const ulong arb_hypgeom_gamma_tab_limbs[];

/* arb's 3456-bit table of gamma (54 limbs on 64-bit systems) */
#define EULER_TAB_LIMBS (3456 / FLINT_BITS)

void
_mp_real_const_euler_compute(mp_real_t res, slong n)
{
    if (n + 1 < EULER_TAB_LIMBS)
    {
        /* gamma in [1/2, 1): the table as a fraction, one ulp */
        _mp_real_set_mpn_2exp(res, arb_hypgeom_gamma_tab_limbs + EULER_TAB_LIMBS,
            EULER_TAB_LIMBS, -FLINT_BITS * EULER_TAB_LIMBS);
        _mp_real_add_error_ulps_at(res, 1.0, -EULER_TAB_LIMBS);
        return;
    }

    _mp_real_const_euler_tune(res, n, -1);
}
