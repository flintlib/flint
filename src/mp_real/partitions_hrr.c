/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    p(n) by the Hardy-Ramanujan-Rademacher formula

        p(n) = sum_{k=1}^N A_k(n) sqrt(3/k) 4/(24n-1) (cosh(C/k) - sinh(C/k)/(C/k)) + R_N,
        C = pi sqrt(24n-1)/6,

    A_k(n) = prefactor sqrt(p/q) prod cos(pi p_i/q_i) from
    arith_hrr_expsum_factored, |R_N| from partitions_rademacher_bound.
    Term k is needed to about 2^-guard absolutely: prec_k = log2 |term
    bound| + guard bits.  The terms of up to about 200 bits are summed in
    batches of dfloat balls (double to quad-double) with the vector
    functions, the others in mp_real balls at ceil(prec_k / 64) limbs
    (when dfloat is not supported, all of them).

    Threads (flint_get_num_threads, from n = 10^8): the constants pi, C,
    exp(C) are computed once with the whole budget (the parallelism is
    inside the library: binary splitting, the bit-burst cascade, the
    FFT multiplications); then the k = 1..N are cut into ranges of about
    equal estimated cost, taken from a queue by the threads, each range
    summed into its own accumulator (shared read-only constants, views
    instead of truncated copies); the accumulators are added in order.
*/

#include <math.h>
#include <float.h>
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "arith.h"
#include "partitions.h"
#include "partitions/impl.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "dfloat.h"
#include "mp_real.h"
#include "impl.h"

#define PI 3.141592653589793238462643
#define INV_LOG2 (1.44269504088896340735992468 + 1e-12)
#define HRR_D (1.2424533248940001551 + 1e-12)
#define MIN_PREC 20
/* the largest exp whose kernel needs no tables (measured: a first call
   at 8 limbs builds tables for 140-450 us, against 1-2 us warm and 3-4 us
   for a root) */
#define EXP_NOTAB_LIMBS 7

/* the bounds of partitions_hrr_sum_arb (the arb implementation this
   replaces): the number of prime factors of k */
static const ulong primorial_tab[] = {
    1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690, 223092870,
#if FLINT_BITS == 64
    UWORD(6469693230), UWORD(200560490130), UWORD(7420738134810),
    UWORD(304250263527210), UWORD(13082761331670030), UWORD(614889782588491410)
#endif
};

static int
bound_primes(ulong k)
{
    unsigned int i;
    for (i = 0; i < sizeof(primorial_tab) / sizeof(ulong); i++)
        if (k <= primorial_tab[i])
            return i;
    return i;
}

/* ceil(log2(x)) as bitcount(floor(x)) (one too many for a power of
   two); but the cast overflows from x = 2^63
   (8 N 26 sqrt(n) for k = 1 reaches it near n = 2.4 10^17: then the
   conversion is undefined, on x86-64 giving 64 instead of 72 bits at
   n = 10^20), so large x by the logarithm, rounded up */
static slong
log2_ceil(double x)
{
    if (x < 0x1p62)
        return FLINT_BIT_COUNT((slong) x);
    return (slong) ceil(log2(x) * (1.0 + 1e-15)) + 1;
}

/* exp(C/k) as exp(C)^(1/k) where exp would build tables: the root is
   table-free and costs less than even a warm table exp from 64 limbs for
   k up to 256 (measured; a cold table costs 0.25-20 ms at 16-512 limbs,
   too much for the few calls of one p(n)) */
static int
use_root(slong k, slong w)
{
    return k == 1 || w > EXP_NOTAB_LIMBS;
}

/* ------------------------------------------------------------------------
   The batched dfloat path.  A term needing prec_k bits goes to the
   smallest format dL (L = 1..4: double to quad-double balls) with
   53 L - log2 z - df_slack[L] >= prec_k; the rest (and everything when
   dfloat is not supported) goes to mp_real.  The slack covers the bits
   the ball evaluation loses: the radius of z = C/k amplified by exp
   (log2 z), and two dozen roundings (df_slack). */

/* bits lost beyond log2 z, measured over 300 random n from 2000 to
   10^10: 3.6, 6.1, 7.6, 8.5; with a 1.5-bit margin */
static const double df_slack[5] = { 0, 5.1, 7.6, 9.1, 10.0 };
/* the smallest batches worth vectorizing, when mp_real is set up anyway */
static slong df_min_batch[5] = { 0, 1, 2, 4, WORD_MAX };
#define DF_LEVELS 4
#ifndef DF_CHUNK
#define DF_CHUNK 256
#endif

#define DF_PQ_MAX 0x1p53     /* products sp sq exact in a double */

typedef struct
{
    slong m, ma;            /* terms */
    double * k, * pre, * sp, * sq;
    slong * f0;
    int * nf;
    slong c, ca;            /* factors */
    double * fa, * fb;
    unsigned char * kind;   /* 0: cos(pi fa/fb), 1: sin(pi fa/fb) */
}
df_batch;

static void
df_batch_init(df_batch * B)
{
    memset(B, 0, sizeof(df_batch));
}

static void
df_batch_clear(df_batch * B)
{
    flint_free(B->k); flint_free(B->pre); flint_free(B->sp); flint_free(B->sq);
    flint_free(B->f0); flint_free(B->nf);
    flint_free(B->fa); flint_free(B->fb); flint_free(B->kind);
}

/* the term k into B; returns 0 (nothing done) if its integers do not fit
   the doubles exactly */
static int
df_push(df_batch * B, const trig_prod_t prod, slong k)
{
    ulong p = prod->sqrt_p * 3, q = prod->sqrt_q * k, g;
    double pre = 4.0 * prod->prefactor;
    slong i;

    /* (the bounds as doubles: they exceed 32-bit words) */
    if ((double) prod->sqrt_p >= 0x1p50 || (double) prod->sqrt_q >= 0x1p50
        || (double) k >= 0x1p50 || q / k != prod->sqrt_q)
        return 0;
    g = n_gcd(p, q);
    p /= g;
    q /= g;
    if ((double) p >= 0x1p53 || (double) q >= 0x1p53)
        return 0;
    for (i = 0; i < prod->n; i++)
        if ((double) prod->cos_q[i] >= 0x1p50)
            return 0;

    if (B->m == B->ma)
    {
        B->ma = FLINT_MAX(16, 2 * B->ma);
        B->k = flint_realloc(B->k, B->ma * sizeof(double));
        B->pre = flint_realloc(B->pre, B->ma * sizeof(double));
        B->sp = flint_realloc(B->sp, B->ma * sizeof(double));
        B->sq = flint_realloc(B->sq, B->ma * sizeof(double));
        B->f0 = flint_realloc(B->f0, B->ma * sizeof(slong));
        B->nf = flint_realloc(B->nf, B->ma * sizeof(int));
    }
    if (B->c + prod->n > B->ca)
    {
        B->ca = FLINT_MAX(B->c + prod->n, FLINT_MAX(32, 2 * B->ca));
        B->fa = flint_realloc(B->fa, B->ca * sizeof(double));
        B->fb = flint_realloc(B->fb, B->ca * sizeof(double));
        B->kind = flint_realloc(B->kind, B->ca * sizeof(unsigned char));
    }

    B->f0[B->m] = B->c;
    B->nf[B->m] = prod->n;
    for (i = 0; i < prod->n; i++)
    {
        /* cos(pi t/qq): |p| mod 2qq, reflections to 0 <= t <= qq/2, then
           a cosine for t/qq <= 1/4 and sin(pi (qq - 2t)/(2qq)) above */
        ulong qq = prod->cos_q[i], t = (ulong) FLINT_ABS(prod->cos_p[i]) % (2 * qq);
        if (t > qq)
            t = 2 * qq - t;
        if (2 * t > qq)
        {
            t = qq - t;
            pre = -pre;
        }
        if (4 * t <= qq)
        {
            B->fa[B->c] = (double) t;
            B->fb[B->c] = (double) qq;
            B->kind[B->c] = 0;
        }
        else
        {
            B->fa[B->c] = (double) (qq - 2 * t);
            B->fb[B->c] = (double) (2 * qq);
            B->kind[B->c] = 1;
        }
        B->c++;
    }
    B->k[B->m] = (double) k;
    B->pre[B->m] = pre;
    B->sp[B->m] = (double) p;
    B->sq[B->m] = (double) q;
    B->m++;
    return 1;
}

/* The sum of a batch's terms as an mp_real ball, in the dfloat balls with L
   components through the generic ring and its vector methods.  For term
   i, with z = C/k:

       term = pre sqrt(sp/sq) prod_j f_j (e^z (1 - 1/z) + e^-z (1 + 1/z)) / (2 (24n - 1)),

   f_j = cos(pi a/b) or sin(pi a/b) with 0 <= a/b <= 1/4 (the reductions
   of df_push: the sines take the factors near pi/2, where a cosine would
   only be known absolutely), pre = 4 x prefactor x the signs of the
   reductions.  Ball arithmetic throughout, so the result encloses the
   exact sum whatever the rounding.  The work vectors of DF_CHUNK terms at
   a time stay in the caches. */

#define E(v, i) GR_ENTRY((v), (i), sz)

/* the constants pi, C = pi sqrt(24n - 1)/6, 1/C, -1/C, 1/(2 (24n - 1)),
   once per p(n) in quad-double balls (K[0..4] in the context c4); each
   format converts them */
static void
df_consts(gr_ptr K, const fmpz_t n24, gr_ctx_t c4)
{
    slong sz = c4->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr pi = K, C = E(K, 1), invC = E(K, 2), ninvC = E(K, 3), inv24h = E(K, 4);

    status |= gr_pi(pi, c4);
    status |= gr_set_fmpz(inv24h, n24, c4);
    status |= gr_sqrt(C, inv24h, c4);
    status |= gr_mul_2exp_si(inv24h, inv24h, 1, c4);
    status |= gr_inv(inv24h, inv24h, c4);
    status |= gr_mul(C, C, pi, c4);
    status |= gr_div_ui(C, C, 6, c4);
    status |= gr_inv(invC, C, c4);
    status |= gr_neg(ninvC, invC, c4);
    GR_MUST_SUCCEED(status);
}

static int
df_eval(mp_real_t res, const df_batch * B, int L, gr_srcptr K4, gr_ctx_t c4, double pqmax)
{
    gr_ctx_t ctx, dctx;
    slong sz, sz4 = c4->sizeof_elem, t0, m, c, cf0, i, j, r, cnt, n0, mmax, cmax, maxnf, nbig;
    gr_ptr K, pi, C, invC, ninvC, inv24h, tot;
    gr_ptr A, z, kb, w1, w2, e, ei, ones, Crep, arg, den, val;
    double * dbuf, * dk, * dpq, * dsp, * dsr, * dpre, * da, * db;
    slong * pos, * perm, * idx, * nfc;
    int status = GR_SUCCESS, ok;

    GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, L, DFLOAT_BALL));
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(dctx, 1, 0));    /* plain doubles */
    sz = ctx->sizeof_elem;

    K = gr_heap_init_vec(7, ctx);
    for (i = 0; i < 5; i++)
        status |= gr_set_other(E(K, i), GR_ENTRY(K4, i, sz4), c4, ctx);
    pi = K; C = E(K, 1); invC = E(K, 2); ninvC = E(K, 3); inv24h = E(K, 4);
    tot = E(K, 6);
    status |= gr_zero(tot, ctx);

    mmax = FLINT_MIN(B->m, DF_CHUNK);
    for (cmax = 0, maxnf = 0, t0 = 0; t0 < B->m; t0 += DF_CHUNK)
    {
        m = FLINT_MIN(DF_CHUNK, B->m - t0);
        cmax = FLINT_MAX(cmax, B->f0[t0 + m - 1] + B->nf[t0 + m - 1] - B->f0[t0]);
    }
    for (i = 0; i < B->m; i++)
        maxnf = FLINT_MAX(maxnf, B->nf[i]);

    A = gr_heap_init_vec(9 * mmax + 3 * cmax + 1, ctx);
    z = E(A, mmax); kb = E(z, mmax); w1 = E(kb, mmax); w2 = E(w1, mmax);
    e = E(w2, mmax); ei = E(e, mmax); ones = E(ei, mmax); Crep = E(ones, mmax);
    arg = E(Crep, mmax); den = E(arg, cmax); val = E(den, cmax);
    dbuf = flint_malloc((5 * mmax + 2 * cmax + 1) * sizeof(double));
    dk = dbuf; dpq = dk + mmax; dsp = dpq + mmax; dpre = dsp + mmax;
    dsr = dpre + mmax; da = dsr + mmax; db = da + cmax;
    pos = flint_malloc((cmax + 3 * mmax + maxnf + 3) * sizeof(slong));
    perm = pos + cmax + 1; idx = perm + mmax; nfc = idx + mmax;
    for (i = 0; i < mmax; i++)
    {
        status |= gr_one(E(ones, i), ctx);
        status |= gr_set(E(Crep, i), C, ctx);
    }

    for (t0 = 0; t0 < B->m; t0 += DF_CHUNK)
    {
        m = FLINT_MIN(DF_CHUNK, B->m - t0);
        cf0 = B->f0[t0];
        c = B->f0[t0 + m - 1] + B->nf[t0 + m - 1] - cf0;

        /* the factors: cosines first, then sines, arguments pi a/b */
        for (n0 = 0, j = 0; j < c; j++)
            n0 += (B->kind[cf0 + j] == 0);
        {
            slong i0 = 0, i1 = n0;
            for (j = 0; j < c; j++)
            {
                slong p = B->kind[cf0 + j] ? i1++ : i0++;
                pos[j] = p;
                da[p] = B->fa[cf0 + j];
                db[p] = B->fb[cf0 + j];
            }
        }
        if (c > 0)
        {
            status |= _gr_vec_set_other(arg, da, dctx, c, ctx);
            status |= _gr_vec_set_other(den, db, dctx, c, ctx);
            status |= _gr_vec_div(arg, arg, den, c, ctx);
            status |= _gr_vec_mul_scalar(arg, arg, c, pi, ctx);
            status |= _gr_vec_cos(val, arg, n0, ctx);
            status |= _gr_vec_sin(E(val, n0), E(arg, n0), c - n0, ctx);
        }

        /* the terms of the chunk ordered by decreasing number of factors
           (the order of a sum is free): the r-th factors are then those
           of a prefix, multiplied in by one gather and one vector
           product per r */
        for (r = 0; r <= maxnf; r++)
            nfc[r] = 0;
        for (i = 0; i < m; i++)
            nfc[B->nf[t0 + i]]++;
        for (cnt = 0, r = maxnf; r >= 0; r--)
        {
            slong t = nfc[r];
            nfc[r] = cnt;
            cnt += t;
        }
        for (i = 0; i < m; i++)
            perm[nfc[B->nf[t0 + i]]++] = t0 + i;

        /* A = pre sqrt(sp/sq) prod f_j, sqrt(sp/sq) = sp rsqrt(sp sq) with
           sp sq an exact double below 2^53: no division, and exact
           inputs keep the SIMD blocks of rsqrt uniform (quotients sq/sp,
           a mix of exact and inexact balls, send whole blocks to scalar).
           Larger products (from k near 3 10^7, i.e. most dfloat terms
           from n near 10^17 on) as sqrt(sp) rsqrt(sq), also on exact
           inputs: the chunk then takes one vector square root more, with
           ones in place of sp for the others (rather than a scalar
           division and square root per term) */
        nbig = 0;
        for (i = 0; i < m; i++)
        {
            slong t = perm[i];
            double pq = B->sp[t] * B->sq[t];
            int small = (pq < pqmax);
            dk[i] = B->k[t];
            dpq[i] = small ? pq : B->sq[t];
            dsp[i] = small ? B->sp[t] : 1.0;
            dsr[i] = small ? 1.0 : B->sp[t];
            dpre[i] = B->pre[t];
            nbig += !small;
        }
        status |= _gr_vec_set_other(w1, dpq, dctx, m, ctx);
        status |= _gr_vec_set_other(w2, dsp, dctx, m, ctx);
        status |= _gr_vec_rsqrt(A, w1, m, ctx);
        status |= _gr_vec_mul(A, A, w2, m, ctx);
        if (nbig)
        {
            status |= _gr_vec_set_other(w1, dsr, dctx, m, ctx);
            status |= _gr_vec_sqrt(w2, w1, m, ctx);
            status |= _gr_vec_mul(A, A, w2, m, ctx);
        }
        status |= _gr_vec_set_other(w1, dpre, dctx, m, ctx);
        status |= _gr_vec_mul(A, A, w1, m, ctx);
        for (r = 0; ; r++)
        {
            for (cnt = 0; cnt < m && B->nf[perm[cnt]] > r; cnt++)
                idx[cnt] = pos[B->f0[perm[cnt]] + r - cf0];
            if (cnt == 0)
                break;
            status |= _gr_vec_gather(w1, val, idx, cnt, ctx);
            status |= _gr_vec_mul(A, A, w1, cnt, ctx);
        }

        /* U = e^z (1 - 1/z) + e^-z (1 + 1/z), z = C/k, 1/z = k/C; e^-z as
           1/e^z (a division costs less than an exponential; its relative
           radius is that of e^z, and the term is e^-2z smaller) */
        status |= _gr_vec_set_other(kb, dk, dctx, m, ctx);
        status |= _gr_vec_div(z, Crep, kb, m, ctx);
        status |= _gr_vec_set(w1, ones, m, ctx);
        status |= _gr_vec_addmul_scalar(w1, kb, m, ninvC, ctx);
        status |= _gr_vec_set(w2, ones, m, ctx);
        status |= _gr_vec_addmul_scalar(w2, kb, m, invC, ctx);
        status |= _gr_vec_exp(e, z, m, ctx);
        status |= _gr_vec_div(ei, ones, e, m, ctx);
        status |= _gr_vec_mul(e, e, w1, m, ctx);
        status |= _gr_vec_mul(ei, ei, w2, m, ctx);
        status |= _gr_vec_add(e, e, ei, m, ctx);

        /* the sum of the terms A U, by a dot product */
        status |= _gr_vec_dot(tot, tot, 0, A, e, m, ctx);
    }
    /* times 1/(2 (24n - 1)) */
    status |= gr_mul(tot, tot, inv24h, ctx);

    /* a failure (not expected: every divisor and square root argument is
       bounded away from zero) is reported to the caller, which redoes the
       sum in mp_real */
    ok = (status == GR_SUCCESS && ((const double *) tot)[L] <= DBL_MAX);
    if (ok)
        mp_real_set_dfloat(res, (const double *) tot, L, ((const double *) tot)[L]);

    gr_heap_clear_vec(A, 9 * mmax + 3 * cmax + 1, ctx);
    flint_free(dbuf);
    flint_free(pos);
    gr_heap_clear_vec(K, 7, ctx);
    gr_ctx_clear(dctx);
    gr_ctx_clear(ctx);
    return ok;
}

#undef E

/* ------------------------------------------------------------------------
   The mp_real constants at W limbs, shared read-only by the tasks: set up
   once (with the whole thread budget for the multiplications and binary
   splittings inside), and only when some term needs mp_real. */

typedef struct
{
    slong W;
    ulong n24ui;        /* 24n - 1 if it fits a word, else 0 (then inv24) */
    mp_real_t n24, C, invC, inv24, E, one;
}
mp_consts;

static void
mp_consts_init(mp_consts * K)
{
    mp_real_init(K->n24); mp_real_init(K->C); mp_real_init(K->invC);
    mp_real_init(K->inv24); mp_real_init(K->E); mp_real_init(K->one);
}

static void
mp_consts_clear(mp_consts * K)
{
    mp_real_clear(K->n24); mp_real_clear(K->C); mp_real_clear(K->invC);
    mp_real_clear(K->inv24); mp_real_clear(K->E); mp_real_clear(K->one);
}

/* E = exp(C) without tables: through log 2 when it is already cached
   (the reduction then costs nothing); else the AGM Newton-Taylor step
   from 16000 limbs on one thread (1.0 - 1.2 times faster than the
   squarings at 18k - 183k limbs) and the squarings below.  With
   threads, the bit-burst cascade inside the squarings runs as parallel
   tasks while the AGM is a chain of full-precision products and square
   roots, so the squarings are faster up to larger sizes (at 2 threads
   0.75 - 0.9 times the AGM time at 58k - 183k limbs, equal at 578k);
   but their tasks hold about twice the memory of a pi computation, so
   the AGM takes over from EXPC_AGM_LIMBS_MT (1 MB numbers) where the
   memory starts to matter. */
#ifndef EXPC_AGM_LIMBS
#define EXPC_AGM_LIMBS 16000
#endif
#ifndef EXPC_AGM_LIMBS_MT
#define EXPC_AGM_LIMBS_MT 131072
#endif

static void
_hrr_exp_C(mp_real_t E, const mp_real_t C, slong w)
{
    if (w + 2 <= MP_REAL_CONST_STATIC_N
        || _mp_real_const_cached_limbs(MP_REAL_CONST_ID_LOG2) >= w + 2)
        mp_real_exp_notab_log2(E, C, w);
    else if (w >= ((flint_get_num_threads() == 1) ? EXPC_AGM_LIMBS : EXPC_AGM_LIMBS_MT))
        mp_real_exp_agm(E, C, w);
    else
        mp_real_exp_notab_squaring(E, C, w);
}

static void
mp_consts_setup(mp_consts * K, const fmpz_t n24f, slong W, int test)
{
    mp_real_t pi4, t;
    mp_real_init(pi4);
    mp_real_init(t);
    K->W = W;
    mp_real_set_ui(K->one, 1);
    /* n24 = 24n - 1, C = (pi/4) sqrt(n24) (2/3), E = exp(C), 1/C, 1/n24;
       C at one more limb: it reaches 2^32 for n near 10^18 */
    mp_real_set_fmpz(K->n24, n24f);
    mp_real_const_pi4(pi4, W + 1, 1);
    mp_real_sqrt(t, K->n24, W + 1);
    mp_real_mul(K->C, pi4, t, W + 1);
    mp_real_clear(pi4);
    mp_real_clear(t);
    mp_real_mul_ui(K->C, K->C, 2, W + 1);
    mp_real_div_ui(K->C, K->C, 3, W + 1);
    _hrr_exp_C(K->E, K->C, W);
    mp_real_div(K->invC, K->one, K->C, W);
    /* the terms divide by 24n - 1 when it is a word (from n near
       7.7 10^17 they multiply by its inverse), and C itself is only read
       by the terms at most EXP_NOTAB_LIMBS long (the others take roots
       of E): neither is held at W limbs through the sum */
    K->n24ui = (fmpz_abs_fits_ui(n24f) && !test) ? fmpz_get_ui(n24f) : 0;
    if (K->n24ui == 0)
        mp_real_div(K->inv24, K->one, K->n24, W);
    mp_real_set_trunc(K->C, K->C, EXP_NOTAB_LIMBS + 4);
}

/* v = a read-only view of x truncated to its top n limbs, sharing the
   limbs (never to be written or cleared): the dropped tail and the old
   radius are each below one ulp of the new bottom limb.  The scalar
   operations read their whole input, so the terms take views of the
   constants at their own precision instead of copies. */
static void
_hrr_view(mp_real_t v, const mp_real_t x, slong n)
{
    *v = *x;
    if (x->size > n)
    {
        v->d = x->d + (x->size - n);
        v->size = n;
        v->err = 1 + (x->err != 0);
    }
}

/* the temporaries of a task */
typedef struct
{
    mp_real_t T, t, e, ei, U;
}
mp_work;

static void
mp_work_init(mp_work * S)
{
    mp_real_init(S->T); mp_real_init(S->t); mp_real_init(S->e);
    mp_real_init(S->ei); mp_real_init(S->U);
}

static void
mp_work_clear(mp_work * S)
{
    mp_real_clear(S->T); mp_real_clear(S->t); mp_real_clear(S->e);
    mp_real_clear(S->ei); mp_real_clear(S->U);
}

/* acc += the term k at bits bits */
static void
mp_term(mp_real_t acc, mp_work * S, const mp_consts * K, const trig_prod_t prod,
    slong k, slong bits, double Cd)
{
    slong wk, i, wi;
    double zd = Cd / k, lg = zd * INV_LOG2;     /* log2 e^z = z / log 2 */
    mp_real_t Ct, invCt, inv24t, Et;
    mp_real_srcptr e;

    /* one limb more than the bits: a ball's top limb may have up to
       63 leading zero bits, and its relative precision is 64 wk less those */
    wk = (bits + FLINT_BITS - 1) / FLINT_BITS + 1;
    _hrr_view(Ct, K->C, wk + 2);
    _hrr_view(invCt, K->invC, wk + 2);
    _hrr_view(inv24t, K->inv24, wk + 2);
    _hrr_view(Et, K->E, wk + 2);

    /* T = 4 prefactor sqrt(3 sqrt_p / (k sqrt_q)) prod cos / n24 */
    {
        ulong p = prod->sqrt_p * 3, q = prod->sqrt_q * k, g, hi, pq;

        g = n_gcd(p, q);
        p /= g;
        q /= g;
        umul_ppmm(hi, pq, p, q);
        if (hi == 0)
        {
            mp_real_rsqrt_ui(S->T, pq, wk);
            mp_real_mul_ui(S->T, S->T, p, wk);
        }
        else
        {
            mp_real_rsqrt_ui(S->T, p, wk);
            mp_real_mul_ui(S->T, S->T, p, wk);
            mp_real_rsqrt_ui(S->t, q, wk);
            mp_real_mul(S->T, S->T, S->t, wk);
        }
        mp_real_mul_ui(S->T, S->T, 4 * (ulong) FLINT_ABS(prod->prefactor), wk);
        if (prod->prefactor < 0)
            mp_real_neg(S->T, S->T);
        for (i = 0; i < prod->n; i++)
        {
            /* relative accuracy of ceil(bits / 64) limbs (wk less the
               guard limb of the products) */
            mp_real_sin_cos_pi_ui_div_ui(NULL, S->t,
                (ulong) FLINT_ABS(prod->cos_p[i]), prod->cos_q[i], FLINT_MAX(wk - 1, 1));
            mp_real_mul(S->T, S->T, S->t, wk);
        }
        if (K->n24ui != 0)
            mp_real_div_ui(S->T, S->T, K->n24ui, wk);
        else
            mp_real_mul(S->T, S->T, inv24t, wk);
    }

    /* U = cosh z - sinh z / z = (e (1 - 1/z) + ei (1 + 1/z)) / 2,
       e = exp(z), ei = 1/e (relatively e^-2z: at the limbs where it is
       visible, or a bound) */
    if (k == 1)
        e = Et;
    else
    {
        if (use_root(k, wk))
            mp_real_root_ui(S->e, Et, k, wk);
        else
        {
            mp_real_div_ui(S->t, Ct, k, wk);
            mp_real_exp_bits(S->e, S->t, bits + 16);
        }
        e = S->e;
    }
    mp_real_mul_ui(S->t, invCt, k, wk);        /* 1/z */
    mp_real_sub(S->U, K->one, S->t, wk);
    mp_real_mul(S->U, S->U, e, wk);
    wi = wk - (slong) ((2 * lg - 2) / FLINT_BITS);
    if (wi <= 0)
    {
        /* ei (1 + 1/z) < 2 e^-z < 2^(1 - floor(lg)) */
        mp_real_add_error_2exp_si(S->U, 1 - (slong) floor(lg));
    }
    else
    {
        mp_real_div(S->ei, K->one, e, wi);
        mp_real_add(S->t, K->one, S->t, wi);
        mp_real_mul(S->ei, S->ei, S->t, wi);
        mp_real_add(S->U, S->U, S->ei, wk);
    }
    mp_real_mul_2exp_si(S->U, S->U, -1);
    mp_real_mul(S->T, S->T, S->U, wk);

    /* an in-place addition of a short term costs O(its length) (the
       accumulator's bottom limb sits at the terms' absolute precision) */
    mp_real_add(acc, acc, S->T, K->W);
}


/* ------------------------------------------------------------------------
   The sum as tasks: ranges of contiguous k, each computing its factored
   sums, its terms (every term in its own format: mp_real, or dfloat
   batches per format, evaluated every DF_FLUSH terms so that the memory
   stays bounded whatever N) and their sum in its own accumulator; the
   accumulators are added in task order at the end.  On one thread the
   whole range is one task.  With threads (and enough work), the range
   is cut into about TASKS_PER_THREAD tasks per thread of about equal
   estimated cost, in increasing k, i.e. by decreasing cost per term (the
   longest mp_real terms first), for the dynamic queue of
   _mp_real_parallel_tasks: the threads take them as they become free,
   which balances the costly early terms against the cheap tail. */

#ifndef DF_FLUSH
#define DF_FLUSH 4096
#endif
#ifndef TASKS_PER_THREAD
#define TASKS_PER_THREAD 16
#endif
/* threads only from this n */
#ifndef PAR_MIN_N
#define PAR_MIN_N 1e6
#endif
/* the estimated costs in ns, for the task sizes: the factored sum per
   k; a dfloat term per format; an mp_real term at w limbs (fitted to
   the task times at n = 10^12, within a factor 1.5 from 5 to 30000
   limbs) */
#define COST_FACT 450.0
static const double cost_df[DF_LEVELS + 1] = { 0, 100, 200, 400, 1000 };

static double
cost_mp(slong w)
{
    return 1000.0 + 90.0 * w * log2(w + 1.0);
}

/* the per-term precision (the bound of the arb implementation) and the
   format, from constants
   of n computed once (one logarithm per term) */
typedef struct
{
    double Cd, l2Cd, lt, s26;
    slong N;
    int df;
}
hrr_bounds;

static void
hrr_bounds_init(hrr_bounds * H, double nd, slong N, int df)
{
    H->Cd = PI * sqrt(24 * nd - 1) / 6;
    H->l2Cd = log2(H->Cd);
    H->lt = HRR_D - log(24.0 * nd - 1);
    H->s26 = 26 * sqrt(nd);
    H->N = N;
    H->df = df;
}

/* bits: the precision of term k (at least MIN_PREC);
   returns the format: the smallest dfloat format whose precision less
   the losses covers the bits, for z = C/k >= 2 (below, the cancellation
   in cosh z - sinh z / z eats into the slack), else 0 (mp_real) */
static int
term_info(slong * bits, const hrr_bounds * H, slong k)
{
    double lk = log2((double) k), z = H->Cd / k, zs;
    slong b;
    int L;

    b = (slong) ((z + H->lt + 0.5 * 0.69314718055994530942 * lk) * INV_LOG2)
        + log2_ceil(8 * H->N * (H->s26 / k + 7 * bound_primes(k) + 22));
    *bits = b = FLINT_MAX(b, MIN_PREC);
    if (!H->df || z < 2.0)
        return 0;
    zs = H->l2Cd - lk;
    for (L = 1; L <= DF_LEVELS && 53 * L - zs - df_slack[L] < b; L++)
        ;
    return (L <= DF_LEVELS) ? L : 0;
}

static slong
term_limbs(slong bits)
{
    /* one limb more than the bits: a ball's top limb may have up to
       63 leading zero bits */
    return (bits + FLINT_BITS - 1) / FLINT_BITS + 1;
}

/* the order of the formats along k: mp_real (5) for the first terms,
   then the dfloat formats from the widest; nonincreasing in k up to
   rounding effects of the bounds (the tasks only use it for their cost
   estimates and the counts below, every term finds its own format) */
static int
term_rank(const hrr_bounds * H, slong k)
{
    slong b;
    int L = term_info(&b, H, k);
    return (L == 0) ? DF_LEVELS + 1 : L;
}

static double
term_cost_mp(const hrr_bounds * H, slong k)
{
    slong b;
    term_info(&b, H, k);
    return COST_FACT + cost_mp(term_limbs(b));
}

typedef struct
{
    slong k0, k1;
}
hrr_task;

typedef struct
{
    const fmpz * n;
    hrr_bounds H;
    double Cd;
    slong W;
    int demote[DF_LEVELS + 1];  /* formats sent to mp_real */
    const mp_consts * K;        /* NULL when not set up */
    gr_srcptr K4;
    gr_ctx_struct * c4;
    double pqmax;               /* DF_PQ_MAX, or 1 to test the other path */
    const hrr_task * task;
    mp_real_struct * acc;
    int * status;               /* per task: HRR_OK, HRR_NEED_MP, HRR_DF_FAIL */
}
hrr_job;

#define HRR_OK 0
#define HRR_NEED_MP 1
#define HRR_DF_FAIL 2

static int
df_flush(mp_real_t acc, df_batch * B, int L, const hrr_job * J, mp_real_t t)
{
    int ok = df_eval(t, B, L, J->K4, J->c4, J->pqmax);
    if (ok)
        mp_real_add(acc, acc, t, J->W);
    B->m = B->c = 0;
    return ok;
}

static void
hrr_worker(slong i, void * arg)
{
    const hrr_job * J = (const hrr_job *) arg;
    const hrr_task * T = J->task + i;
    mp_real_struct * acc = J->acc + i;
    mp_work S;
    df_batch B[DF_LEVELS + 1];
    trig_prod_t prod;
    mp_real_t t;
    slong k, bits, wk, ws = 0;
    int L;

    mp_work_init(&S);
    for (L = 1; L <= DF_LEVELS; L++)
        df_batch_init(B + L);
    mp_real_init(t);
    mp_real_zero(acc);
    J->status[i] = HRR_OK;

    for (k = T->k0; k < T->k1; k++)
    {
        trig_prod_init(prod);
        arith_hrr_expsum_factored(prod, k, fmpz_fdiv_ui(J->n, k));
        if (prod->prefactor == 0)
            continue;
        L = term_info(&bits, &J->H, k);
        if (J->demote[L])
            L = 0;

        if (L > 0 && df_push(B + L, prod, k))
        {
            if (B[L].m >= DF_FLUSH && !df_flush(acc, B + L, L, J, t))
            {
                J->status[i] = HRR_DF_FAIL;
                break;
            }
            continue;
        }
        /* an mp_real term (or a dfloat one whose integers do not fit the
           doubles: not seen) without the constants: the task is redone
           once they are there */
        if (J->K == NULL)
        {
            J->status[i] = HRR_NEED_MP;
            break;
        }
        /* the temporaries of a much longer earlier term (k = 1 at W
           limbs) are released rather than held through the task */
        wk = term_limbs(bits);
        if (ws > 2 * wk + 8)
        {
            mp_work_clear(&S);
            mp_work_init(&S);
            ws = 0;
        }
        ws = FLINT_MAX(ws, wk);
        mp_term(acc, &S, J->K, prod, k, bits, J->Cd);
    }
    for (L = 1; L <= DF_LEVELS; L++)
    {
        if (B[L].m && J->status[i] == HRR_OK && !df_flush(acc, B + L, L, J, t))
            J->status[i] = HRR_DF_FAIL;
        df_batch_clear(B + L);
    }

    mp_real_clear(t);
    mp_work_clear(&S);
}

/* appends the tasks for [k0, k1) of rank r, of about target cost each */
static void
add_tasks(hrr_task ** tasks, slong * ntasks, slong * alloc, slong k0, slong k1,
    int r, double target, const hrr_bounds * H)
{
    slong k = k0;
    while (k < k1)
    {
        slong a = k;
        if (r <= DF_LEVELS)
        {
            /* equal pieces */
            slong len = (slong) FLINT_MAX(1.0, target / (COST_FACT + cost_df[r]));
            k = FLINT_MIN(k1, a + len);
        }
        else
        {
            double c = 0.0;
            while (k < k1 && (c < target || k == a))
            {
                c += term_cost_mp(H, k);
                k++;
            }
        }
        if (*ntasks == *alloc)
        {
            *alloc = FLINT_MAX(16, 2 * *alloc);
            *tasks = flint_realloc(*tasks, *alloc * sizeof(hrr_task));
        }
        (*tasks)[*ntasks].k0 = a;
        (*tasks)[*ntasks].k1 = k;
        (*ntasks)++;
    }
}

static int
hrr_attempt(mp_real_t res, const fmpz_t n, int df, int test)
{
    double nd, Cd;
    slong N, k, bits1, W, i, ntasks, talloc, nt, workers;
    slong cnt[DF_LEVELS + 1], kr[DF_LEVELS + 3];
    int L, r, need;
    fmpz_t n24f;
    mp_consts K;
    int have_K = 0;
    hrr_task * tasks = NULL;
    hrr_job J;
    gr_ctx_t c4;
    gr_ptr K4 = NULL;
    arf_t bound;
    mp_real_t X;
    int ok;

    nd = fmpz_get_d(n);
    N = partitions_hrr_needed_terms(nd);
    Cd = PI * sqrt(24 * nd - 1) / 6;
    nt = flint_get_num_threads();

    fmpz_init(n24f);
    fmpz_mul_ui(n24f, n, 24);
    fmpz_sub_ui(n24f, n24f, 1);

    hrr_bounds_init(&J.H, nd, N, df);
    term_info(&bits1, &J.H, 1);
    W = term_limbs(FLINT_MAX(bits1, 64));

    /* the ranges of the formats: k in [kr[r + 1], kr[r]) has rank r,
       by bisection on the (nonincreasing) rank */
    kr[DF_LEVELS + 2] = 1;
    kr[0] = N + 1;
    for (r = DF_LEVELS + 1; r >= 1; r--)
    {
        /* the first k with rank < r, at or after the range above */
        slong lo = kr[r + 1], hi = N + 1;
        while (lo < hi)
        {
            slong mid = lo + (hi - lo) / 2;
            if (term_rank(&J.H, mid) < r)
                hi = mid;
            else
                lo = mid + 1;
        }
        kr[r] = lo;
    }
    for (L = 0; L <= DF_LEVELS; L++)
        cnt[L] = kr[(L == 0) ? DF_LEVELS + 1 : L] - kr[(L == 0) ? DF_LEVELS + 1 + 1 : L + 1];

    /* Small batches of the wider formats cost more than mp_real terms
       (measured per term against an mp_real term at the same precision:
       d2b from 2 terms, d3b from 4; d4b never, its scalar operations are
       2-3 times slower than the mp_real ones), so they go to mp_real,
       unless no term needs mp_real at all: then the whole mp_real setup
       (pi, exp(C)) is saved. */
    J.demote[0] = 0;
    need = (cnt[0] > 0);
    for (L = 1; L <= DF_LEVELS; L++)
    {
        J.demote[L] = (cnt[0] > 0 && cnt[L] < df_min_batch[L]);
        need |= (J.demote[L] && cnt[L] > 0);
    }

    /* the tasks */
    ntasks = talloc = 0;
    workers = (nt >= 2 && (nd >= PAR_MIN_N || test)) ? nt : 1;
    if (workers == 1 && !test)
    {
        add_tasks(&tasks, &ntasks, &talloc, 1, N + 1, 1, 2.0 * N * (COST_FACT + cost_df[1]), &J.H);
    }
    else
    {
        double total = 0.0, target;
        for (r = DF_LEVELS + 1; r >= 1; r--)
        {
            if (r <= DF_LEVELS)
                total += (kr[r] - kr[r + 1]) * (COST_FACT + cost_df[r]);
            else
                for (k = kr[r + 1]; k < kr[r]; k++)
                    total += term_cost_mp(&J.H, k);
        }
        /* (tested on one thread too, as several tasks) */
        target = total / (TASKS_PER_THREAD * FLINT_MAX(workers, 2));
        for (r = DF_LEVELS + 1; r >= 1; r--)
            add_tasks(&tasks, &ntasks, &talloc, kr[r + 1], kr[r],
                (r <= DF_LEVELS && J.demote[r]) ? DF_LEVELS + 1 : r, target, &J.H);
    }

    /* the constants: mp_real with the whole thread budget, dfloat once */
    mp_consts_init(&K);
    if (need)
    {
        mp_consts_setup(&K, n24f, W, test);
        have_K = 1;
    }
    if (df)
    {
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(c4, 4, DFLOAT_BALL));
        K4 = gr_heap_init_vec(5, c4);
        df_consts(K4, n24f, c4);
    }

    J.n = n; J.Cd = Cd; J.W = W;
    J.pqmax = test ? 1.0 : DF_PQ_MAX;
    J.K = have_K ? &K : NULL;
    J.K4 = K4; J.c4 = df ? c4 : NULL;
    J.task = tasks;
    J.acc = flint_malloc(ntasks * sizeof(mp_real_struct));
    J.status = flint_malloc(ntasks * sizeof(int));
    for (i = 0; i < ntasks; i++)
        mp_real_init(J.acc + i);

    _mp_real_parallel_tasks_max(hrr_worker, &J, ntasks, workers);

    /* the tasks that met an mp_real term without the constants (not
       seen), redone with them; a dfloat failure (not seen either) fails
       the attempt, which the caller redoes in mp_real */
    ok = 1;
    for (i = 0; i < ntasks; i++)
    {
        if (J.status[i] == HRR_NEED_MP)
        {
            if (!have_K)
            {
                mp_consts_setup(&K, n24f, W, test);
                have_K = 1;
                J.K = &K;
            }
            hrr_worker(i, &J);
        }
        if (J.status[i] != HRR_OK)
            ok = 0;
    }

    /* the constants are no longer needed (at huge precision: two
       full-length numbers less through the conversions below) */
    mp_consts_clear(&K);
    mp_consts_init(&K);

    /* the sum of the accumulators, in task order */
    mp_real_init(X);
    for (i = 0; i < ntasks; i++)
    {
        mp_real_add(X, X, J.acc + i, W);
        mp_real_clear(J.acc + i);
    }

    /* the remainder bound, and the integer, kept as an mp_real (an fmpz
       would overflow the int limb count of mpz beyond n near 1.4e21):
       the bound (an arf) goes in as the sum of the powers of two of its
       30-bit upper bound, each exact above the sum's ulp */
    arf_init(bound);
    partitions_rademacher_bound(bound, n, N);
    {
        mag_t mb;
        slong e, i;
        mag_init(mb);
        arf_get_mag(mb, bound);
        if (!mag_is_zero(mb))
        {
            e = fmpz_get_si(MAG_EXPREF(mb)) - MAG_BITS;
            for (i = 0; i < MAG_BITS; i++)
                if ((MAG_MAN(mb) >> i) & 1)
                    mp_real_add_error_2exp_si(X, e + i);
        }
        mag_clear(mb);
    }
    /* in place: no second full-length number */
    ok = ok && mp_real_unique_integer(X, X);
    if (ok)
        mp_real_swap(res, X);
    mp_real_clear(X);

    arf_clear(bound);
    flint_free(J.acc);
    flint_free(J.status);
    flint_free(tasks);
    if (df)
    {
        gr_heap_clear_vec(K4, 5, c4);
        gr_ctx_clear(c4);
    }
    mp_consts_clear(&K);
    fmpz_clear(n24f);
    return ok;
}

void
_mp_real_partitions_hrr(mp_real_t res, const fmpz_t n, int flags)
{
    int df = !(flags & MP_REAL_PARTITIONS_NO_DFLOAT) && dfloat_is_supported();
    int test = (flags & MP_REAL_PARTITIONS_TEST) != 0;

    /* (the bounds need n >= 2) */
    if (fmpz_cmp_ui(n, 2) < 0)
    {
        if (fmpz_sgn(n) >= 0)
            mp_real_set_ui(res, 1);
        else
            mp_real_zero(res);
        return;
    }

    if (hrr_attempt(res, n, df, test))
        return;
    /* not expected: the dfloat radii are bounded by the slack (and a
       dfloat failure is not expected either); redo everything in mp_real,
       whose precisions carry a whole guard limb */
    if (df && hrr_attempt(res, n, 0, test))
        return;
    flint_throw(FLINT_ERROR, "mp_real_partitions_hrr: no unique integer\n");
}

/* the lookup table, then the pentagonal recurrence (exact below 417 on
   64-bit: p(416) < 2^64 <= p(417)); it costs less than the HRR sum up
   to about 500 (0.05 us at 128, 6.8 us at 417, 9.5 us at 500, against
   7.4 - 9.4 us for the sum, which has a floor of about 7 us), so no
   intermediate method pays off */
#define PARTITIONS_REC_MAX ((FLINT_BITS == 64) ? 417 : 128)

void
mp_real_partitions_hrr(mp_real_t res, ulong nhi, ulong nlo)
{
    if (nhi == 0 && nlo < 128)
    {
        mp_real_set_ui(res, partitions_lookup[nlo]);
    }
    else if (nhi == 0 && nlo < PARTITIONS_REC_MAX)
    {
        nn_ptr v;
        TMP_INIT;
        TMP_START;
        v = TMP_ALLOC((nlo + 1) * sizeof(ulong));
        _partitions_vec_ui(v, nlo + 1);
        mp_real_set_ui(res, v[nlo]);
        TMP_END;
    }
    else
    {
        fmpz_t n;
        fmpz_init(n);
        fmpz_set_uiui(n, nhi, nlo);
        _mp_real_partitions_hrr(res, n, 0);
        fmpz_clear(n);
    }
}
