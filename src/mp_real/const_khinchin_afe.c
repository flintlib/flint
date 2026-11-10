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
#include "fmpz.h"
#include "fmpz_mat.h"
#include "mag.h"
#include "thread_support.h"
#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif
#include "mp_real.h"
#include "impl.h"

/* Khinchin's constant by a theta-function splitting (an approximate
   functional equation), in p^(3/2 + o(1)) bit operations.

   With the dilogarithm series of Bailey, Borwein and Crandall,

       log K log 2 = log^2 2 + Li2(-1/2) + (1/2) sum_(k>=2) (-1)^k Li2(4/k^2),

   and Li2(4/K) = int_0^inf e^(-K u) G(u) du (K > 4), G(u) = E(4u)/u,
   E(x) = sum_(n>=1) x^n / (n n!), the alternating sum over k >= 3 is an
   integral of a theta function against G; splitting it at u0 = 4^-a
   and transforming theta_4 on (0, u0] gives, for every u0 > 0,

       log K log 2 = log^2 2 + Li2(-1/2) + pi^2/12 + T/2,
       T = sum_(k>=3) (-1)^k I_k + sum_(m>=0) J_m + M,

       I_k = int_(u0)^inf e^(-k^2 u) G(u) du,
       J_m = sqrt(pi) int_0^(u0) u^(-1/2) e^(-pi^2 (m+1/2)^2/u) G(u) du,
       M   = int_0^(u0) (e^-u - e^(-4u) - 1/2) G(u) du.

   For p bits, about sqrt(p log 2 / u0) k-terms and sqrt(p u0 log 2)/pi
   dual terms are needed.  A k-term is a rational series times
   e^(-k^2 u0), summed by binary splitting (khinchin_afe_kterm.c); a dual
   term is the value of a holonomic function of c = pi^2 (m + 1/2)^2,
   evaluated by the bit-burst algorithm (khinchin_afe_dual.c); M is a
   rational series.  Each term costs O(M(p) log^2 p).

   The terms are independent; they run as parallel tasks, the dual terms
   first (the most expensive ones), at most maxdual of them at once
   since a dual term takes about four times the memory of a k-term
   (mostly transforms, see linrec.c), the other threads
   taking chunks of consecutive k-terms meanwhile (consecutive, so that
   e^(-k^2 u0) can be propagated by multiplications).  The finished
   terms are added into one fixed-point integer with an integer error
   count, so that the memory is that of the running terms and the result
   does not depend on the order in which they finish.

   The working precisions of the terms (relative to their sizes) are
   lowered by about log2(e) times their exponential decay; this needs
   log2(e) to more digits than the margins: the 5-digit 1.4427 is off by
   5e-6 times x bits for a term of size e^-x, x up to p log 2 for p bits,
   which exceeds the 32 guard bits of a k-term from about 3 10^6 digits. */

/* the main term *************************************************************/

/* M = M_1 - M_4 - M_0/2, M_a = int_0^u0 e^(-a u) G(u) du; with
   w^(a)(u) = e^(-au) E(4u) = sum w_N u^N,
   N w_N = -a w_(N-1) + ((4-a)^N - (-a)^N)/N!, so M_a = sum_(N>=1) w_N u0^N / N
   (a > 0) and M_0 = sum 4^N u0^N / (N^2 N!).  The state at index n is
   u0^n (w1_n, w4_n, 3^n/n!, (-1)^n/n!, (-4)^n/n!, 4^n/n!), and the sum
   takes (w1_n - w4_n)/n - e4_n/(2 n^2) of the new state; over the
   denominator 4 n^5 u0d. */
static void
_main_leaf(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong n, const void * vctx)
{
    ulong u0d = *(const ulong *) vctx;
    /* the transition M (over 2 n^3) as small multiples of n and n^2 */
    static const struct { int i, j, c, deg; } ent[] = {
        { 0, 0, -2, 2 }, { 0, 2, 6, 1 }, { 0, 3, 2, 1 },
        { 1, 1, -8, 2 }, { 1, 4, 8, 1 },
        { 2, 2, 6, 2 }, { 3, 3, -2, 2 }, { 4, 4, -8, 2 }, { 5, 5, 8, 2 } };
    fmpz_t m, nn, n2;
    slong k;

    fmpz_init(m);
    fmpz_init_set_ui(nn, n);
    fmpz_init(n2);
    fmpz_mul_ui(n2, nn, 2 * (ulong) n);          /* 2 n^2 */

    for (k = 0; k < (slong) (sizeof(ent) / sizeof(ent[0])); k++)
    {
        /* the state rows, times 2 n^2 */
        fmpz_set_si(m, ent[k].c);
        fmpz_mul_ui(m, m, n);
        if (ent[k].deg == 2)
            fmpz_mul_ui(m, m, n);
        fmpz_mul(fmpz_mat_entry(T, ent[k].i, ent[k].j), m, n2);
        /* the sum row: (2n, -2n, 0, 0, 0, -1) M */
        if (ent[k].i == 0)
            fmpz_addmul_ui(fmpz_mat_entry(T, 6, ent[k].j), m, 2 * (ulong) n);
        else if (ent[k].i == 1)
            fmpz_submul_ui(fmpz_mat_entry(T, 6, ent[k].j), m, 2 * (ulong) n);
        else if (ent[k].i == 5)
            fmpz_sub(fmpz_mat_entry(T, 6, ent[k].j), fmpz_mat_entry(T, 6, ent[k].j), m);
    }

    /* q = 2 n^3 2 n^2 u0d */
    fmpz_mul(q, n2, n2);
    fmpz_mul_ui(q, q, n);
    fmpz_mul_ui(q, q, u0d);

    fmpz_clear(m);
    fmpz_clear(nn);
    fmpz_clear(n2);
}

static void
_main_term(mp_real_t res, slong a, ulong u0d, slong prec)
{
    _mp_real_linrec_struct rec;
    mp_real_struct x0[6];
    mag_t t, f;
    double lf = 0.0, u0 = ldexp(1.0, -2 * a);
    slong N = 10, i, n = mp_real_prec_bits(prec + 32);

    /* the coefficients of e^(-au) E(4u) are bounded by those of e^((a+4)u),
       so the terms are at most 5 (8 u0)^N / N!, and their sum from N on at
       most 10 (8 u0)^N / N! when 8 u0 / (N + 1) <= 1/2 */
    while (lf < prec + 40 || N + 1 < 16 * u0)
    {
        N++;
        lf += log2((double) N) - 3 + 2 * a;
    }

    rec.d = 6;
    rec.s = 1;
    rec.g = 0;
    rec.leaf = _main_leaf;
    rec.bottom = NULL;
    rec.bottom_len = 0;
    rec.ctx = &u0d;
    for (i = 0; i < 6; i++)
    {
        mp_real_init(x0 + i);
        mp_real_set_ui(x0 + i, i >= 2);
    }
    mp_real_zero(res);
    _mp_real_linrec_eval(x0, res, NULL, &rec, 1, N, NULL, 1, n + 1);

    mag_init(t);
    mag_init(f);
    mag_set_ui_2exp_si(t, 1, 3 - 2 * a);
    mag_pow_ui(t, t, N);
    mag_rfac_ui(f, N);
    mag_mul(t, t, f);
    mag_mul_ui(t, t, 10);
    _mp_real_add_error_mag(res, t);
    mag_clear(t);
    mag_clear(f);
    for (i = 0; i < 6; i++)
        mp_real_clear(x0 + i);
}

/* log^2 2 + Li2(-1/2) + pi^2/12, with Li2(-1/2) = -Li2(1/3) - log^2(3/2) / 2
   and Li2(1/3) = Li2(1/9) / 6 + pi^2/18 - log^2 3 / 6 */
static void
_constant_part(mp_real_t res, slong n)
{
    /* Li2(1/9) = sum_(k>=1) P(k)/Q(k) prod_(j<k) R(j)/Q(j), Q = 9 k^2,
       R = k^2, P = 1 */
    static const int64_t P[1] = { 1 }, Q[3] = { 0, 0, 9 }, R[3] = { 0, 0, 1 };
    mp_real_t li, pi2, l2, l3, t;

    mp_real_init(li); mp_real_init(pi2); mp_real_init(l2); mp_real_init(l3); mp_real_init(t);

    mp_real_hypgeom_series_int64(li, 1, 1, 0, 1, P, 1, Q, 3, R, 3, n);
    mp_real_const_pi4(pi2, n, 1);
    mp_real_mul(pi2, pi2, pi2, n);
    mp_real_mul_2exp_si(pi2, pi2, 4);                   /* pi^2 */
    mp_real_const_log2(l2, n, 1);
    mp_real_set_ui(l3, 3);
    mp_real_log_bits(l3, l3, FLINT_BITS * n);

    /* Li2(1/3) */
    mp_real_div_ui(li, li, 6, n);
    mp_real_div_ui(t, pi2, 18, n);
    mp_real_add(li, li, t, n);
    mp_real_mul(t, l3, l3, n);
    mp_real_div_ui(t, t, 6, n);
    mp_real_sub(li, li, t, n);
    /* Li2(-1/2) = -Li2(1/3) - log^2(3/2)/2 */
    mp_real_sub(t, l3, l2, n);
    mp_real_mul(t, t, t, n);
    mp_real_mul_2exp_si(t, t, -1);
    mp_real_add(li, li, t, n);
    mp_real_neg(res, li);
    mp_real_mul(t, l2, l2, n);
    mp_real_add(res, res, t, n);
    mp_real_div_ui(t, pi2, 12, n);
    mp_real_add(res, res, t, n);

    mp_real_clear(li); mp_real_clear(pi2); mp_real_clear(l2); mp_real_clear(l3); mp_real_clear(t);
}

/* scheduling ****************************************************************/

typedef struct
{
    const _mp_real_khinchin_afe_struct * ctx;
    slong mmax, nk, kc;          /* dual terms, k-terms (k >= 3), k-chunks */
    /* the sum of the terms as acc 2^acc_e within accerr 2^acc_e: exact
       integer additions, so that the result does not depend on the
       order in which the terms finish, and no term is kept */
    fmpz_t acc, accerr;
    slong acc_e;
    slong next_d, next_k, dual_active, maxdual;
    int maxdual_hard;            /* maxdual given: a cap to the end */
    slong done_d, done_k;
    double t0, t0_d, t0_k, t_last;   /* t_last: the last report */
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
    pthread_cond_t cond;
#endif
}
_sched_struct;

static void
_lock(_sched_struct * S)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&S->mutex);
#endif
}

static void
_unlock(_sched_struct * S)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&S->mutex);
#endif
}

static void
_report(_sched_struct * S, int dual, slong idx, double t)
{
    double now = _mp_real_verbose_time(), el, eta;
    slong done, tot;

    if (!mp_real_get_verbose())
        return;
    if (dual)
    {
        done = S->done_d;
        tot = S->mmax;
        el = now - S->t0_d;
    }
    else
    {
        done = S->done_k;
        tot = S->nk;
        el = now - S->t0_k;
    }
    eta = (done < tot) ? el / done * (tot - done) : 0.0;
    /* at level 1, at most one line per interval, throttled apart from
       the progress lines of the evaluations inside the terms (which
       would otherwise crowd these out) */
    if (mp_real_get_verbose() >= 2)
        _mp_real_log("khinchin afe: %s %wd: %.3f s; %wd/%wd %s, phase %.0f s, ETA %.0f s",
            dual ? "dual term m =" : "k-terms to k =", idx, t, done, tot,
            dual ? "dual terms" : "k-terms", el, eta);
    else if (now - S->t_last >= MP_REAL_VERBOSE_INTERVAL || done == tot)
    {
        _mp_real_log("khinchin afe: %wd/%wd dual terms, %wd/%wd k-terms; %s ETA %.0f s",
            S->done_d, S->mmax, S->done_k, S->nk, dual ? "dual" : "k", eta);
        S->t_last = now;
    }
}

/* x as an integer r 2^e within u 2^e (u counts the radius and the
   floor) */
static void
_to_units(fmpz_t r, fmpz_t u, const mp_real_t x, slong e)
{
    mag_t m;

    mag_init(m);
    _mp_real_get_fmpz_fixed(r, m, x, e);
    mag_mul_2exp_si(m, m, -e);
    mag_get_fmpz(u, m);
    mag_clear(m);
}

static void
_dual_unit(_sched_struct * S, slong m)
{
    const _mp_real_khinchin_afe_struct * ctx = S->ctx;
    slong cn = mp_real_prec_bits(2 * ctx->prec + 128);
    mp_real_t c, J;
    fmpz_t r, u;
    double t = _mp_real_verbose_time();

    /* c = pi^2 (m + 1/2)^2 */
    mp_real_init(c);
    mp_real_init(J);
    fmpz_init(r);
    fmpz_init(u);
    mp_real_const_pi4(c, cn, 1);
    mp_real_mul_ui(c, c, 2 * (ulong) m + 1, cn);
    mp_real_mul_2exp_si(c, c, 1);
    mp_real_mul(c, c, c, cn);
    _mp_real_khinchin_afe_dual(J, c, ctx);
    _to_units(r, u, J, S->acc_e);
    mp_real_clear(c);
    mp_real_clear(J);

    _lock(S);
    fmpz_add(S->acc, S->acc, r);
    fmpz_add(S->accerr, S->accerr, u);
    S->done_d++;
    S->dual_active--;
    _report(S, 1, m, _mp_real_verbose_time() - t);
#if FLINT_USES_PTHREAD
    pthread_cond_broadcast(&S->cond);
#endif
    _unlock(S);
    fmpz_clear(r);
    fmpz_clear(u);
}

static void
_k_unit(_sched_struct * S, slong j)
{
    const _mp_real_khinchin_afe_struct * ctx = S->ctx;
    slong j0 = j * S->nk / S->kc, j1 = (j + 1) * S->nk / S->kc, i;
    mp_real_t ek, fk, g, x, sum;
    fmpz_t r, u;
    double t = _mp_real_verbose_time();

    mp_real_init(ek);
    mp_real_init(fk);
    mp_real_init(g);
    mp_real_init(x);
    mp_real_init(sum);
    fmpz_init(r);
    fmpz_init(u);

    /* e^(-(k+1)^2 u0) = e^(-k^2 u0) e^(-(2k+1) u0),
       e^(-(2k+3) u0) = e^(-(2k+1) u0) e^(-2 u0) */
    for (i = j0; i < j1; i++)
    {
        ulong k = 3 + i, K = k * k;
        slong wpk = ctx->prec - (slong) ((K - 4.0) * ldexp(1.0, -2 * ctx->a)
            * 1.4426950408889634) + 128;

        wpk = FLINT_MAX(wpk, 128);
        if (i == j0)
        {
            /* by powering e^-u0: about 3 log2(k) multiplications in
               all, under half the cost of one exponential at 3 10^5
               digits; the relative error grows by about the exponent
               K (under 2^30 up to 10^7 digits), within the 192 extra
               bits of e^-u0 */
            slong pn = mp_real_prec_bits(wpk + 64);
            mp_real_pow_ui(ek, ctx->eu, K, pn);
            mp_real_pow_ui(fk, ctx->eu, 2 * k + 1, pn);
            mp_real_mul(g, ctx->eu, ctx->eu, pn);
        }
        else
        {
            mp_real_mul(ek, ek, fk, mp_real_prec_bits(wpk));
            mp_real_mul(fk, fk, g, mp_real_prec_bits(wpk));
        }
        _mp_real_khinchin_afe_kterm(x, k, ctx, ek);
        if (k & 1)
            mp_real_sub(sum, sum, x, mp_real_prec_bits(ctx->prec));
        else
            mp_real_add(sum, sum, x, mp_real_prec_bits(ctx->prec));
    }

    _to_units(r, u, sum, S->acc_e);
    mp_real_clear(ek);
    mp_real_clear(fk);
    mp_real_clear(g);
    mp_real_clear(x);
    mp_real_clear(sum);

    _lock(S);
    fmpz_add(S->acc, S->acc, r);
    fmpz_add(S->accerr, S->accerr, u);
    S->done_k += j1 - j0;
    _report(S, 0, 2 + j1, _mp_real_verbose_time() - t);
    _unlock(S);
    fmpz_clear(r);
    fmpz_clear(u);
}

/* a task takes the next dual term while fewer than maxdual run (or, with
   the default maxdual, when no k-chunk is left), else the next chunk of
   k-terms; with only dual terms left and all slots taken, it waits for
   one to finish */
static void
_unit_task(slong i, void * arg)
{
    _sched_struct * S = (_sched_struct *) arg;
    slong m = -1, j = -1;

    (void) i;
    _lock(S);
    for (;;)
    {
        /* the default cap only holds while k-chunks remain to be
           taken, so that the threads do not idle at the end */
        if (S->next_d < S->mmax && (S->dual_active < S->maxdual
                || (!S->maxdual_hard && S->next_k >= S->kc)))
        {
            m = S->next_d++;
            S->dual_active++;
            if (S->t0_d < 0)
                S->t0_d = _mp_real_verbose_time();
            break;
        }
        if (S->next_k < S->kc)
        {
            j = S->next_k++;
            if (S->t0_k < 0)
                S->t0_k = _mp_real_verbose_time();
            break;
        }
        if (S->next_d >= S->mmax)
            break;
#if FLINT_USES_PTHREAD
        pthread_cond_wait(&S->cond, &S->mutex);
#else
        break;
#endif
    }
    _unlock(S);

    if (m >= 0)
        _dual_unit(S, m);
    else if (j >= 0)
        _k_unit(S, j);
}

/* the driver ****************************************************************/

#define KHINCHIN_AFE_DEFAULT_A 2

/* T = log K = S / log 2 with the terms summed to absolute accuracy
   about 2^-prec */
static void
_khinchin_afe_log(mp_real_t T, slong n, slong a, slong maxdual, slong prec)
{
    _mp_real_khinchin_afe_struct ctx;
    _sched_struct S;
    mp_real_t t;
    mag_t tl;
    slong wn, kmax, m, nt;
    double u0;

    wn = mp_real_prec_bits(prec + 32);
    u0 = ldexp(1.0, -2 * a);

    ctx.a = a;
    ctx.u0d = UWORD(1) << (2 * a);
    ctx.prec = prec;
    mp_real_init(ctx.g1);
    mp_real_init(ctx.e4);
    mp_real_init(ctx.eu);
    mp_real_init(t);
    mag_init(tl);

    /* g1 = E(4 u0) = sum_(k>=1) (4 u0)^k / (k k!): P = 1, Q = d k^2, R = k
       with 4 u0 = 1/d, d = 4^(a - 1) */
    {
        int64_t P[1] = { 1 }, Q[3] = { 0, 0, 0 }, R[2] = { 0, 1 };
        Q[2] = (int64_t) 1 << (2 * (a - 1));
        mp_real_hypgeom_series_int64(ctx.g1, 1, 1, 0, 1, P, 1, Q, 3, R, 2, wn + 1);
    }
    mp_real_set_ui(ctx.e4, 1);
    mp_real_mul_2exp_si(ctx.e4, ctx.e4, 2 - 2 * a);
    mp_real_exp_bits(ctx.e4, ctx.e4, FLINT_BITS * (wn + 1));
    mp_real_set_ui(ctx.eu, 1);
    mp_real_mul_2exp_si(ctx.eu, ctx.eu, -2 * a);
    mp_real_neg(ctx.eu, ctx.eu);
    mp_real_exp_bits(ctx.eu, ctx.eu, prec + 192);

    /* the terms */
    kmax = (slong) sqrt(prec * 0.6931471805599453 / u0 + 4) + 2;
    /* the bound for the omitted k-terms below needs (2 k + 1) u0 >= 0.7
       for k > kmax (only matters for small u0 at low precision) */
    while ((2 * (ulong) kmax + 3) * 10 < 7 * ctx.u0d)
        kmax++;
    for (m = 0; ; m++)
        if (pow(3.14159265358979 * (m + 0.5), 2) / u0 * 1.4426950408889634 > prec + 10)
            break;

    if ((ulong) kmax >= (UWORD(1) << (FLINT_BITS / 2)) - 1)
        flint_throw(FLINT_ERROR, "_mp_real_const_khinchin_afe: precision too high\n");

    nt = flint_get_num_threads();
    S.ctx = &ctx;
    S.mmax = m;
    S.nk = kmax - 2;
    S.kc = FLINT_MIN(S.nk, FLINT_MAX(64, 16 * nt));
    fmpz_init(S.acc);
    fmpz_init(S.accerr);
    S.acc_e = -(prec + 64);
    S.next_d = S.next_k = S.dual_active = 0;
    S.done_d = S.done_k = 0;
    S.maxdual = (maxdual > 0) ? maxdual : FLINT_MAX(1, (nt + 1) / 2);
    S.maxdual_hard = (maxdual > 0);
    S.t0 = _mp_real_verbose_time();
    S.t0_d = S.t0_k = -1.0;
    S.t_last = -1e300;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&S.mutex, NULL);
    pthread_cond_init(&S.cond, NULL);
#endif

    if (mp_real_get_verbose())
        _mp_real_log("khinchin afe: %wd limbs, u0 = 4^-%wd, %wd k-terms (k <= %wd) in %wd chunks, %wd dual terms (at most %wd at once)",
            n, a, S.nk, kmax, S.kc, S.mmax, S.maxdual);

    _mp_real_parallel_tasks(_unit_task, &S, S.mmax + S.kc);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&S.mutex);
    pthread_cond_destroy(&S.cond);
#endif

    /* T = the sum of the terms */
    mp_real_set_fmpz(T, S.acc);
    mp_real_mul_2exp_si(T, T, S.acc_e);
    mag_set_fmpz(tl, S.accerr);
    mag_mul_2exp_si(tl, tl, S.acc_e);
    _mp_real_add_error_mag(T, tl);
    fmpz_clear(S.acc);
    fmpz_clear(S.accerr);

    /* the omitted terms.  k > kmax: 0 < I_k <= e^(-(K-4) u0) / ((K-4) u0),
       the ratios at most e^(-(2k+1) u0) <= 1/2, so the sum is at most twice
       the first.  m >= mmax: 0 < J_m <= sqrt(pi) G(u0) e^(-c/u0)
       int_0^u0 u^(-1/2) du = 2 sqrt(pi u0) G(u0) e^(-c/u0), G(u0) = g1/u0,
       also at most twice the first. */
    {
        ulong k1 = kmax + 1, K1 = k1 * k1;
        mag_t x, y, p;

        mag_init(x);
        mag_init(y);
        mag_init(p);
        mag_set_ui_lower(x, K1 - 4);
        mag_mul_2exp_si(x, x, -2 * a);
        mag_expinv(tl, x);
        mag_div_ui(tl, tl, K1 - 4);
        mag_mul_2exp_si(tl, tl, 2 * a + 1);

        mag_const_pi_lower(x);
        mag_mul_ui_lower(x, x, 2 * S.mmax + 1);
        mag_mul_2exp_si(x, x, -1);
        mag_mul_lower(x, x, x);
        mag_mul_2exp_si(x, x, 2 * a);
        mag_expinv(y, x);
        mag_const_pi(p);
        mag_sqrt(p, p);
        mag_mul(y, y, p);
        _mp_real_get_mag(p, ctx.g1);
        mag_mul(y, y, p);
        mag_mul_2exp_si(y, y, a + 2);               /* 4 sqrt(pi) u0^(-1/2) g1 */
        mag_add(tl, tl, y);
        _mp_real_add_error_mag(T, tl);
        mag_clear(x);
        mag_clear(y);
        mag_clear(p);
    }

    /* + M, then S = T/2 + the constant part, K = exp(S / log 2) */
    _main_term(t, a, ctx.u0d, prec);
    mp_real_add(T, T, t, wn);
    mp_real_mul_2exp_si(T, T, -1);
    _constant_part(t, wn);
    mp_real_add(T, T, t, wn);
    mp_real_const_log2(t, wn, 1);
    mp_real_div(T, T, t, wn);

    mp_real_clear(ctx.g1);
    mp_real_clear(ctx.e4);
    mp_real_clear(ctx.eu);
    mp_real_clear(t);
    mag_clear(tl);
}

void
_mp_real_const_khinchin_afe(mp_real_t res, slong n, slong a, slong maxdual)
{
    mp_real_t T;
    mp_real_phase_t ph;

    if (n < 1)
        flint_throw(FLINT_ERROR, "_mp_real_const_khinchin_afe: n < 1\n");
    if (a <= 0)
        a = KHINCHIN_AFE_DEFAULT_A;
    if (a > 7)
        flint_throw(FLINT_ERROR, "_mp_real_const_khinchin_afe: u0 = 4^-a with a > 7\n");

    MP_REAL_PHASE_START(ph, "khinchin afe", n);
    mp_real_init(T);

    /* log K to absolute accuracy 2^-(FLINT_BITS n + 32) gives K to n
       limbs; the terms are summed with 64 more bits */
    _khinchin_afe_log(T, n, a, maxdual, FLINT_BITS * (n + 1) + 32);
    mp_real_exp_bits(res, T, FLINT_BITS * (n + 1));

    mp_real_clear(T);
    MP_REAL_PHASE_END(ph);
}
