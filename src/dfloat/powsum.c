/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fp_contract.h"
#include <math.h>
#include <string.h>
#include "ulong_extras.h"
#include "double_extras.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "dfloat.h"
#include "thread_support.h"

/* The sieved power sum sum_{n <= N} n^(-s) = sum_{n <= N} f(n) over
   generic rings, with runtime dispatch through gr: the terms live in
   a ball ring T (real numbers; complex values are pairs of real
   vectors), the phases t log p in a ball ring P of higher precision,
   and, optionally, the composites and the block sums in a plain
   (non-ball) ring F with an a priori error analysis.  All the bulk
   work goes through vector kernels (_gr_vec_log, _gr_vec_sin_cos,
   _gr_vec_exp, _gr_vec_mul, _gr_vec_add, _gr_vec_set_other,
   _gr_vec_gather, _gr_vec_scatter, ...).

   The algorithm (as acb_dirichlet_powsum_sieved): f is completely
   multiplicative, so only the primes need transcendental functions;
   f(n) = f(d) f(n/d) for the smallest prime factor d of a composite n.
   The odd numbers are processed in blocks of POWSUM_BLOCK (or all of
   them at once when there are fewer, so that small N do not pay for
   the full-size scratch vectors), with a
   segmented sieve giving the smallest prime factor and its cofactor;
   the values f(k) for odd k <= N/3 are kept in a table (interleaved
   real and imaginary parts, so that a random access touches one cache
   line).  The two-phase version (the default for N >= 2048) keeps
   the table only up to M ~ N^(2/3): beyond M, only the roots c (the
   primes and the composites whose cofactor n/d is <= M) are computed,
   each times G(N/c, d), the sum of f(m) over the m <= N/c whose odd
   part has no prime factor above d = spf(c) (a small table; see
   dfloat.rst), which accounts for all n = m c.  The second phase
   runs in parallel, in chunks of a fixed length.  The powers of 2 are handled by the decomposition
   sum_{n <= N} f(n) = sum_j f(2)^j U_j, U_j = sum_{k odd <= N/2^j} f(k),
   evaluated by Horner's rule from the bucket sums u_j = sum of f(k)
   over 2^-(j+1) N < k <= 2^-j N.

   The prime terms: log p in P, the phase reduced modulo 1 as tau log p
   minus a nearby integer (exact), then everything in T:
   p^(-s) = exp(-sigma log p) (cos(2 pi phi) - i sin(2 pi phi)).

   The a priori error bound of the plain variant (F != NULL), with
   g(k) = k^(-sigma_lo) >= |k^(-s)| for all s in the input ball, F(k)
   the computed value and rho(k) = |F(k) - k^(-s)| / g(k):

   * Primes: rho(p) <= (rad Re + rad Im) / (a lower bound of p^(-sigma))
     from the balls; eps = the maximum over all primes.
   * A plain complex product z = a b in F has |z - ab| <= mu |a| |b|
     with mu = sqrt(2) (dm + da + dm da), where dm and da are relative
     error bounds of the plain product and sum in F (relative to
     |x| |y| and |x| + |y|).  By induction, rho(n) <= (1 + eps)^Omega
     (1 + mu)^(Omega - 1) - 1 with Omega = Omega(n) <= log_3 N.
   * The pairwise sum of a piece of at most 2^D terms has an error of
     at most ((1 + da)^D - 1) sum |F(k)| in each of the real and
     imaginary parts, and sum |F(k)| <= (1 + rho) sum g(k).
   * The piece sums are accumulated into the bucket sums, and those
     combined by Horner's rule, in ball arithmetic.  The radius added
     to bucket j (both parts) is (rho + ((1 + da)^D - 1)(1 + rho)) W_j
     with W_j >= sum g(k) over the odd k of the bucket, bounded by
     the number of terms times max(g(a), g(b)) for the buckets with
     few terms and by max(g(a), g(b)) + (1/2) int_a^b g otherwise. */

#ifndef POWSUM_LOG_BLOCK
#define POWSUM_LOG_BLOCK 12
#endif
#define POWSUM_BLOCK (WORD(1) << POWSUM_LOG_BLOCK)

#define ENTRY(v, i, ctx) GR_ENTRY(v, i, (ctx)->sizeof_elem)

typedef struct
{
    gr_ctx_struct * T, * P;
    gr_ptr xp, lp, ph, nv;      /* P */
    gr_ptr fr, lt, c, sn, m;    /* T */
    int critical;               /* sigma = 1/2: p^(-1/2) by _gr_vec_rsqrt */
}
powsum_scratch_struct;

/* the terms p^(-s) for the primes p[0..len) */
static int
powsum_prime_terms(gr_ptr re, gr_ptr im, const ulong * p, slong len,
    gr_srcptr tau, gr_srcptr negsigma, gr_srcptr twopi, powsum_scratch_struct * S)
{
    gr_ctx_struct * T = S->T, * P = S->P;
    int status = GR_SUCCESS, round;
    slong i;

    for (i = 0; i < len; i++)
        status |= gr_set_ui(ENTRY(S->xp, i, P), p[i], P);
    status |= _gr_vec_log(S->lp, S->xp, len, P);
    status |= _gr_vec_mul_scalar(S->ph, S->lp, len, tau, P);
    /* the fractional part: subtract an integer near the midpoint
       (exact).  When the double approximation of a midpoint is itself
       a large integer (|x| >= 2^52), the residue can still be large
       (the rounding of gr_get_d; for expansions, the next component),
       so the reduction is repeated until all residues are small */
    for (round = 0; round < 8; round++)
    {
        int big = 0;
        for (i = 0; i < len; i++)
        {
            double d;
            status |= gr_get_d(&d, ENTRY(S->ph, i, P), P);
            big |= (fabs(d) >= 0x1p52);
            status |= gr_set_d(ENTRY(S->nv, i, P), rint(d), P);
        }
        status |= _gr_vec_sub(S->ph, S->ph, S->nv, len, P);
        if (!big)
            break;
    }
    status |= _gr_vec_set_other(S->fr, S->ph, P, len, T);
    status |= _gr_vec_set_other(S->lt, S->lp, P, len, T);
    status |= _gr_vec_mul_scalar(S->fr, S->fr, len, twopi, T);
    status |= _gr_vec_sin_cos(S->sn, S->c, S->fr, len, T);
    if (S->critical)
    {
        /* p^(-1/2) from the exact p, more accurate than exp(-(log p)/2) */
        status |= _gr_vec_set_other(S->lt, S->xp, P, len, T);
        status |= _gr_vec_rsqrt(S->m, S->lt, len, T);
    }
    else
    {
        status |= _gr_vec_mul_scalar(S->lt, S->lt, len, negsigma, T);
        status |= _gr_vec_exp(S->m, S->lt, len, T);
    }
    status |= _gr_vec_mul(re, S->m, S->c, len, T);
    status |= _gr_vec_mul(im, S->m, S->sn, len, T);
    status |= _gr_vec_neg(im, im, len, T);
    return status;
}

/* (zr, zi) = (ar, ai) (br, bi), scalars; t1, t2 scratch */
static int
powsum_cmul(gr_ptr zr, gr_ptr zi, gr_srcptr ar, gr_srcptr ai, gr_srcptr br, gr_srcptr bi,
    gr_ptr t1, gr_ptr t2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_mul(t1, ar, br, ctx);
    status |= gr_mul(t2, ai, bi, ctx);
    status |= gr_sub(zr, t1, t2, ctx);
    status |= gr_mul(t1, ar, bi, ctx);
    status |= gr_mul(t2, ai, br, ctx);
    status |= gr_add(zi, t1, t2, ctx);
    return status;
}

/* the sum of x[0..len) (destroyed) by pairwise vector additions */
static int
powsum_tree_sum(gr_ptr res, gr_ptr x, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    if (len == 0)
        return gr_zero(res, ctx);
    while (len > 1)
    {
        slong half = len / 2;
        status |= _gr_vec_add(x, x, ENTRY(x, half, ctx), half, ctx);
        if (len & 1)
            status |= gr_set(ENTRY(x, half, ctx), ENTRY(x, len - 1, ctx), ctx);
        len = half + (len & 1);
    }
    status |= gr_set(res, x, ctx);
    return status;
}

/* an upper bound for sum_{k odd, a <= k <= b} k^(-sigma) (a, b odd),
   sigma >= sigma_lo: the number of terms times the largest term, which
   is a^(-sigma_lo) or b^(-sigma_lo). On the dyadic pieces (b <= 2a)
   this is within a factor 2^|sigma_lo| of the sum (20% for sigma =
   1/2), but cheap (mag arithmetic): for the pieces with few terms. */
static double
powsum_weight_bound_count(ulong a, ulong b, const arf_t sigma_lo)
{
    mag_t x, y;
    arf_t r;
    double res;

    if (a > b)
        return 0.0;
    mag_init(x);
    mag_init(y);
    arf_init(r);
    if (arf_sgn(sigma_lo) >= 0)
    {
        /* a^(-sigma) <= exp(-sigma_lo log(a)) */
        arf_get_mag_lower(x, sigma_lo);
        mag_set_ui_lower(y, a);
        mag_log_lower(y, y);
        mag_mul_lower(x, x, y);
        mag_expinv(x, x);
    }
    else
    {
        /* b^(-sigma) <= exp(|sigma_lo| log(b)) */
        arf_get_mag(x, sigma_lo);
        mag_log_ui(y, b);
        mag_mul(x, x, y);
        mag_exp(x, x);
    }
    mag_mul_ui(x, x, (b - a) / 2 + 1);
    arf_set_mag(r, x);
    res = arf_get_d(r, ARF_RND_UP);
    mag_clear(x);
    mag_clear(y);
    arf_clear(r);
    return res;
}

/* the same, sharper for many terms: max(g(a), g(b)) + (1/2) int_a^b g,
   with g(x) = x^(-sigma) */
static double
powsum_weight_bound_integral(ulong a, ulong b, const arb_t sigma)
{
    arb_t A, B, t, u, e;
    arf_t r;
    double res;
    slong prec = 64;

    if (a > b)
        return 0.0;
    arb_init(A); arb_init(B); arb_init(t); arb_init(u); arb_init(e);
    arf_init(r);
    arb_set_ui(A, a);
    arb_set_ui(B, b);
    arb_neg(e, sigma);
    /* max(g(a), g(b)) */
    arb_pow(t, A, e, prec);
    arb_pow(u, B, e, prec);
    arb_union(t, t, u, prec);
    /* (1/2) int_a^b x^(-sigma) dx */
    if (arb_is_one(sigma))
    {
        arb_div(u, B, A, prec);
        arb_log(u, u, prec);
    }
    else
    {
        arb_add_ui(e, e, 1, prec);            /* 1 - sigma */
        arb_pow(u, B, e, prec);
        {
            arb_t v;
            arb_init(v);
            arb_pow(v, A, e, prec);
            arb_sub(u, u, v, prec);
            arb_clear(v);
        }
        arb_div(u, u, e, prec);
    }
    arb_mul_2exp_si(u, u, -1);
    arb_add(t, t, u, prec);
    arb_get_ubound_arf(r, t, prec);
    res = arf_get_d(r, ARF_RND_UP);
    arb_clear(A); arb_clear(B); arb_clear(t); arb_clear(u); arb_clear(e);
    arf_clear(r);
    return res;
}

static double
powsum_weight_bound(ulong a, ulong b, const arb_t sigma_lo)
{
    if (b - a < 256)
        return powsum_weight_bound_count(a, b, arb_midref(sigma_lo));
    else
        return powsum_weight_bound_integral(a, b, sigma_lo);
}

/* an upper bound for k^(-sigma), sigma >= sigma_lo, added to res */
static void
powsum_g_add(mag_t res, ulong k, const arf_t sigma_lo)
{
    mag_t x, y;
    mag_init(x);
    mag_init(y);
    if (arf_sgn(sigma_lo) >= 0)
    {
        arf_get_mag_lower(x, sigma_lo);
        mag_set_ui_lower(y, k);
        mag_log_lower(y, y);
        mag_mul_lower(x, x, y);
        mag_expinv(x, x);
    }
    else
    {
        arf_get_mag(x, sigma_lo);
        mag_log_ui(y, k);
        mag_mul(x, x, y);
        mag_exp(x, x);
    }
    mag_add(res, res, x);
    mag_clear(x);
    mag_clear(y);
}

static double
powsum_mag_get_d_up(const mag_t x)
{
    arf_t r;
    double d;
    arf_init(r);
    arf_set_mag(r, x);
    d = arf_get_d(r, ARF_RND_UP);
    arf_clear(r);
    return d;
}

/* The table of G(z, p) = sum of f(m) over the m = 2^a m' <= z with m'
   odd and (if p != 0) its prime factors at most p, for 1 <= z <= Z,
   as balls in T: the row of G(z) = G(z, 0) (z = 0..Z, at 0..Z), and
   for each odd prime p < Z the row of G(z, p) for p < z <= Z (for
   z <= p, G(z, p) = G(z)). row[p] is the offset of the row of p. */
static int
powsum_G_table(gr_ptr Gre, gr_ptr Gim, slong * row, ulong Z,
    gr_ptr Xre, gr_ptr Xim, ulong * list, slong blk,
    gr_srcptr tau, gr_srcptr negsigma, gr_srcptr twopi, powsum_scratch_struct * S)
{
    gr_ctx_struct * T = S->T;
    ulong * lpf, z, p, k, off;
    gr_ptr fre, fim;
    int status = GR_SUCCESS;

    lpf = flint_malloc(sizeof(ulong) * (Z + 1));
    fre = gr_heap_init_vec(Z + 1, T);
    fim = gr_heap_init_vec(Z + 1, T);

    /* the largest prime factor of the odd part */
    for (z = 0; z <= Z; z++)
        lpf[z] = 1;
    for (p = 3; p <= Z; p += 2)
        if (lpf[p] == 1)
            for (k = p; k <= Z; k += 2 * p)
                lpf[k] = p;
    for (z = 2; z <= Z; z += 2)
        lpf[z] = lpf[z >> flint_ctz(z)];

    /* f(m), m = 1..Z, directly */
    for (z = 1; z <= Z && status == GR_SUCCESS; )
    {
        slong len = FLINT_MIN((slong) (Z - z + 1), blk), i;
        for (i = 0; i < len; i++)
            list[i] = z + i;
        status |= powsum_prime_terms(Xre, Xim, list, len, tau, negsigma, twopi, S);
        status |= _gr_vec_set(ENTRY(fre, z, T), Xre, len, T);
        status |= _gr_vec_set(ENTRY(fim, z, T), Xim, len, T);
        z += len;
    }

    /* the full row */
    status |= gr_zero(Gre, T);
    status |= gr_zero(Gim, T);
    for (z = 1; z <= Z; z++)
    {
        status |= gr_add(ENTRY(Gre, z, T), ENTRY(Gre, z - 1, T), ENTRY(fre, z, T), T);
        status |= gr_add(ENTRY(Gim, z, T), ENTRY(Gim, z - 1, T), ENTRY(fim, z, T), T);
    }

    /* the rows of the odd primes p < Z */
    off = Z + 1;
    for (p = 3; p < Z && status == GR_SUCCESS; p += 2)
    {
        gr_srcptr pr, pi;
        if (lpf[p] != p)
            continue;
        row[p] = off;
        pr = ENTRY(Gre, p, T);
        pi = ENTRY(Gim, p, T);
        for (z = p + 1; z <= Z; z++)
        {
            if (lpf[z] <= p)
            {
                status |= gr_add(ENTRY(Gre, off, T), pr, ENTRY(fre, z, T), T);
                status |= gr_add(ENTRY(Gim, off, T), pi, ENTRY(fim, z, T), T);
            }
            else
            {
                status |= gr_set(ENTRY(Gre, off, T), pr, T);
                status |= gr_set(ENTRY(Gim, off, T), pi, T);
            }
            pr = ENTRY(Gre, off, T);
            pi = ENTRY(Gim, off, T);
            off++;
        }
    }

    gr_heap_clear_vec(fre, Z + 1, T);
    gr_heap_clear_vec(fim, Z + 1, T);
    flint_free(lpf);
    return status;
}

/* the number of entries of the table above */
static slong
powsum_G_length(ulong Z)
{
    slong len = Z + 1;
    ulong p, k;
    char * comp = flint_calloc(Z + 1, 1);
    for (p = 3; p < Z; p += 2)
        if (!comp[p])
        {
            len += Z - p;
            for (k = p * p; k <= Z; k += 2 * p)
                comp[k] = 1;
        }
    flint_free(comp);
    return len;
}

/* ceil(N^(2/3)), at least sqrt(N) + 1 */
static ulong
powsum_default_table_bound(ulong N)
{
    ulong M = (ulong) ceil(pow((double) N, 2.0 / 3.0));
    return FLINT_MAX(M, n_sqrt(N) + 1);
}

/* the per-thread scratch: vectors of the block length, the sieve
   arrays of the segment length, the sieve state */
typedef struct
{
    powsum_scratch_struct S;
    gr_ctx_struct * F;
    gr_ptr Fre, Fim, Ar, Ai, Br, Bi, T1, T2, T3, T4;    /* F */
    gr_ptr Pre, Pim;                                    /* T */
    ulong * pl, * spf, * cof, * snext, * scof;
    slong * ppos, * cpos, * ia, * ib, * even, * gix;
    slong blk, seg, nsmall;
}
powsum_work_struct;

static void
powsum_work_init(powsum_work_struct * W, gr_ctx_t T, gr_ctx_t P, gr_ctx_t F,
    slong blk, slong seg, slong nsmall, int critical)
{
    slong i;
    W->F = F;
    W->blk = blk;
    W->seg = seg;
    W->nsmall = nsmall;
    W->Fre = gr_heap_init_vec(blk, F);
    W->Fim = gr_heap_init_vec(blk, F);
    W->Ar = gr_heap_init_vec(blk, F);
    W->Ai = gr_heap_init_vec(blk, F);
    W->Br = gr_heap_init_vec(blk, F);
    W->Bi = gr_heap_init_vec(blk, F);
    W->T1 = gr_heap_init_vec(blk, F);
    W->T2 = gr_heap_init_vec(blk, F);
    W->T3 = gr_heap_init_vec(blk, F);
    W->T4 = gr_heap_init_vec(blk, F);
    W->Pre = gr_heap_init_vec(blk, T);
    W->Pim = gr_heap_init_vec(blk, T);
    W->S.T = T;
    W->S.P = P;
    W->S.critical = critical;
    W->S.xp = gr_heap_init_vec(blk, P);
    W->S.lp = gr_heap_init_vec(blk, P);
    W->S.ph = gr_heap_init_vec(blk, P);
    W->S.nv = gr_heap_init_vec(blk, P);
    W->S.fr = gr_heap_init_vec(blk, T);
    W->S.lt = gr_heap_init_vec(blk, T);
    W->S.c = gr_heap_init_vec(blk, T);
    W->S.sn = gr_heap_init_vec(blk, T);
    W->S.m = gr_heap_init_vec(blk, T);
    W->pl = flint_malloc(sizeof(ulong) * blk);
    W->spf = flint_malloc(sizeof(ulong) * seg);
    W->cof = flint_malloc(sizeof(ulong) * seg);
    W->snext = flint_malloc(sizeof(ulong) * (nsmall + 1));
    W->scof = flint_malloc(sizeof(ulong) * (nsmall + 1));
    W->ppos = flint_malloc(sizeof(slong) * blk);
    W->cpos = flint_malloc(sizeof(slong) * blk);
    W->ia = flint_malloc(sizeof(slong) * blk);
    W->ib = flint_malloc(sizeof(slong) * blk);
    W->even = flint_malloc(sizeof(slong) * blk);
    W->gix = flint_malloc(sizeof(slong) * blk);
    for (i = 0; i < blk; i++)
        W->even[i] = 2 * i;
}

static void
powsum_work_clear(powsum_work_struct * W)
{
    gr_ctx_struct * F = W->F, * T = W->S.T, * P = W->S.P;
    slong blk = W->blk;
    gr_heap_clear_vec(W->Fre, blk, F);
    gr_heap_clear_vec(W->Fim, blk, F);
    gr_heap_clear_vec(W->Ar, blk, F);
    gr_heap_clear_vec(W->Ai, blk, F);
    gr_heap_clear_vec(W->Br, blk, F);
    gr_heap_clear_vec(W->Bi, blk, F);
    gr_heap_clear_vec(W->T1, blk, F);
    gr_heap_clear_vec(W->T2, blk, F);
    gr_heap_clear_vec(W->T3, blk, F);
    gr_heap_clear_vec(W->T4, blk, F);
    gr_heap_clear_vec(W->Pre, blk, T);
    gr_heap_clear_vec(W->Pim, blk, T);
    gr_heap_clear_vec(W->S.xp, blk, P);
    gr_heap_clear_vec(W->S.lp, blk, P);
    gr_heap_clear_vec(W->S.ph, blk, P);
    gr_heap_clear_vec(W->S.nv, blk, P);
    gr_heap_clear_vec(W->S.fr, blk, T);
    gr_heap_clear_vec(W->S.lt, blk, T);
    gr_heap_clear_vec(W->S.c, blk, T);
    gr_heap_clear_vec(W->S.sn, blk, T);
    gr_heap_clear_vec(W->S.m, blk, T);
    flint_free(W->pl); flint_free(W->spf); flint_free(W->cof);
    flint_free(W->snext); flint_free(W->scof);
    flint_free(W->ppos); flint_free(W->cpos); flint_free(W->ia); flint_free(W->ib);
    flint_free(W->even); flint_free(W->gix);
}

/* the data shared by the ranges */
typedef struct
{
    gr_ctx_struct * T, * P, * F;
    int plain, critical;
    ulong N, Mt, L1, Z;
    slong blk, seg, nsmall;
    const ulong * small;
    gr_ptr tab;                 /* F; written by the first phase only */
    gr_srcptr Gre, Gim;         /* F (the midpoints, for a plain F) */
    const slong * row;
    gr_srcptr tau, negsigma, twopi;
    gr_ptr ure, uim;            /* T: the bucket sums of the first phase */
    /* the second phase, in chunks of clen odd numbers from c0 */
    ulong c0, clen, cN;
    gr_ptr Rre, Rim;            /* T, one per chunk */
    double * eps;               /* one per chunk */
    int * status;               /* one per chunk */
}
powsum_args_struct;

/* The odd numbers in [lo0, hi0] (odd, within the first phase [1, L1]
   or the second (L1, N]): the first phase adds to the bucket sums and
   writes the table; the second adds its sum over the roots to
   (Rre, Rim).  The largest relative error of the prime terms (plain
   F) goes to *eps. The sieve runs over segments of seg odd numbers
   (so that the primes up to sqrt(N) are visited about once per
   sqrt(N) numbers), the vector work over blocks of blk. */
static int
powsum_range(const powsum_args_struct * A, ulong lo0, ulong hi0,
    gr_ptr Rre, gr_ptr Rim, double * eps_out)
{
    gr_ctx_struct * T = A->T, * P = A->P, * F = A->F;
    powsum_work_struct W;
    powsum_scratch_struct * S;
    gr_ptr tab = A->tab;
    gr_ptr t1, t2, s1, s2;
    ulong N = A->N, Mt = A->Mt, L1 = A->L1, slo, shi, lo, hi, k;
    slong j, i, np, nc, jtop = 0, blk = A->blk, seg = A->seg, nsmall = A->nsmall;
    const ulong * small = A->small;
    double eps = 0.0;
    int plain = A->plain, phase2 = (lo0 > L1);
    int status = GR_SUCCESS;

    powsum_work_init(&W, T, P, F, blk, seg, nsmall, A->critical);
    S = &W.S;
    GR_TMP_INIT2(t1, t2, T);
    GR_TMP_INIT2(s1, s2, F);

    /* the sieve state at lo0: for each small prime q, its first odd
       multiple >= max(q^2, lo0) and the cofactor */
    for (j = 0; j < nsmall; j++)
    {
        ulong q = small[j], m;
        if (q * q >= lo0)
            m = q * q;
        else
        {
            m = ((lo0 + q - 1) / q) * q;
            if (m % 2 == 0)
                m += q;
        }
        W.snext[j] = m;
        W.scof[j] = m / q;
    }

    for (slo = lo0; slo <= hi0 && status == GR_SUCCESS; slo = shi + 2)
    {
        slong scnt;

        shi = FLINT_MIN(slo + 2 * (seg - 1), hi0);
        scnt = (shi - slo) / 2 + 1;

        /* the smallest prime factor and its cofactor: the primes are
           marked in decreasing order, so that the smallest one writes
           last (no branch) */
        for (i = 0; i < scnt; i++)
            W.spf[i] = 0;
        while (jtop < nsmall && small[jtop] * small[jtop] <= shi)
            jtop++;
        for (j = jtop - 1; j >= 0; j--)
        {
            ulong q = small[j], m = W.snext[j], c = W.scof[j];
            for (; m <= shi; m += 2 * q, c += 2)
            {
                W.spf[(m - slo) / 2] = q;
                W.cof[(m - slo) / 2] = c;
            }
            W.snext[j] = m;
            W.scof[j] = c;
        }

        for (lo = slo; lo <= shi && status == GR_SUCCESS; lo = hi + 2)
        {
            slong cnt;
            const ulong * spf, * cof;

            hi = FLINT_MIN(lo + 2 * (blk - 1), shi);
            cnt = (hi - lo) / 2 + 1;
            spf = W.spf + (lo - slo) / 2;
            cof = W.cof + (lo - slo) / 2;

            np = nc = 0;
            for (i = 0; i < cnt; i++)
            {
                int isp = (spf[i] == 0);
                W.ppos[np] = i;
                W.pl[np] = lo + 2 * i;
                W.cpos[nc] = i;
                np += isp;
                nc += !isp;
            }
            if (lo == 1)
            {
                /* 1 is not a prime */
                status |= gr_one(W.Fre, F);
                status |= gr_zero(W.Fim, F);
                for (i = 0; i + 1 < np; i++)
                {
                    W.ppos[i] = W.ppos[i + 1];
                    W.pl[i] = W.pl[i + 1];
                }
                np--;
            }

            /* the primes, in ball arithmetic */
            status |= powsum_prime_terms(W.Pre, W.Pim, W.pl, np, A->tau, A->negsigma, A->twopi, S);
            if (plain)
            {
                /* their relative errors (rad Re + rad Im) / (m - rad m),
                   as balls (exact inputs) in T; the tiny radii of the
                   quotients are covered by the slack */
                status |= _gr_vec_get_interval_mid_rad(S->fr, S->lt, W.Pre, np, T);
                status |= _gr_vec_get_interval_mid_rad(S->fr, S->c, W.Pim, np, T);
                status |= _gr_vec_add(S->lt, S->lt, S->c, np, T);
                status |= _gr_vec_get_interval_mid_rad(S->fr, S->c, S->m, np, T);
                status |= _gr_vec_sub(S->fr, S->fr, S->c, np, T);
                if (_gr_vec_div(S->sn, S->lt, S->fr, np, T) != GR_SUCCESS)
                    eps = D_INF;
                for (i = 0; i < np; i++)
                {
                    double r;
                    status |= gr_get_d(&r, ENTRY(S->sn, i, T), T);
                    r = r * (1.0 + 0x1p-40);
                    if (!(r >= 0.0) || !(r <= 1e-3))
                        r = D_INF;
                    eps = FLINT_MAX(eps, r);
                }
            }

            if (phase2)
            {
                /* the roots: the primes and the composites with
                   cofactor <= Mt, in Fre, Fim (primes first), each
                   times G(N/c, d), d = spf(c) (0 for a prime) */
                slong nr, nrc = 0;
                ulong c, z;

                if (plain)
                {
                    status |= _gr_vec_set_other(W.Fre, W.Pre, T, np, F);
                    status |= _gr_vec_set_other(W.Fim, W.Pim, T, np, F);
                }
                else
                {
                    status |= _gr_vec_set(W.Fre, W.Pre, np, F);
                    status |= _gr_vec_set(W.Fim, W.Pim, np, F);
                }
                for (i = 0; i < np; i++)
                {
                    c = W.pl[i];
                    z = (ulong) ((double) N / (double) c);
                    if (z * c > N) z--;
                    else if ((z + 1) * c <= N) z++;
                    W.gix[i] = z;
                }
                for (i = 0; i < nc; i++)
                {
                    ulong d = spf[W.cpos[i]], e = cof[W.cpos[i]];
                    if (e <= Mt)
                    {
                        c = lo + 2 * W.cpos[i];
                        W.ia[nrc] = d - 1;
                        W.ib[nrc] = e - 1;
                        z = (ulong) ((double) N / (double) c);
                        if (z * c > N) z--;
                        else if ((z + 1) * c <= N) z++;
                        W.gix[np + nrc] = (d < z) ? A->row[d] + (slong) (z - d - 1) : (slong) z;
                        nrc++;
                    }
                }
                status |= _gr_vec_gather(W.Ar, tab, W.ia, nrc, F);
                status |= _gr_vec_gather(W.Ai, ENTRY(tab, 1, F), W.ia, nrc, F);
                status |= _gr_vec_gather(W.Br, tab, W.ib, nrc, F);
                status |= _gr_vec_gather(W.Bi, ENTRY(tab, 1, F), W.ib, nrc, F);
                status |= _gr_vec_mul(W.T1, W.Ar, W.Br, nrc, F);
                status |= _gr_vec_mul(W.T2, W.Ai, W.Bi, nrc, F);
                status |= _gr_vec_mul(W.T3, W.Ar, W.Bi, nrc, F);
                status |= _gr_vec_mul(W.T4, W.Ai, W.Br, nrc, F);
                status |= _gr_vec_sub(ENTRY(W.Fre, np, F), W.T1, W.T2, nrc, F);
                status |= _gr_vec_add(ENTRY(W.Fim, np, F), W.T3, W.T4, nrc, F);
                nr = np + nrc;

                status |= _gr_vec_gather(W.Br, A->Gre, W.gix, nr, F);
                status |= _gr_vec_gather(W.Bi, A->Gim, W.gix, nr, F);
                status |= _gr_vec_mul(W.T1, W.Fre, W.Br, nr, F);
                status |= _gr_vec_mul(W.T2, W.Fim, W.Bi, nr, F);
                status |= _gr_vec_mul(W.T3, W.Fre, W.Bi, nr, F);
                status |= _gr_vec_mul(W.T4, W.Fim, W.Br, nr, F);
                status |= _gr_vec_sub(W.T1, W.T1, W.T2, nr, F);
                status |= _gr_vec_add(W.T3, W.T3, W.T4, nr, F);
                status |= powsum_tree_sum(s1, W.T1, nr, F);
                status |= powsum_tree_sum(s2, W.T3, nr, F);
                status |= gr_set_other(t1, s1, F, T);
                status |= gr_set_other(t2, s2, F, T);
                status |= gr_add(Rre, Rre, t1, T);
                status |= gr_add(Rim, Rim, t2, T);
                continue;
            }

            if (plain)
            {
                status |= _gr_vec_set_other(W.T1, W.Pre, T, np, F);
                status |= _gr_vec_scatter(W.Fre, W.ppos, W.T1, np, F);
                status |= _gr_vec_set_other(W.T1, W.Pim, T, np, F);
                status |= _gr_vec_scatter(W.Fim, W.ppos, W.T1, np, F);
            }
            else
            {
                status |= _gr_vec_scatter(W.Fre, W.ppos, W.Pre, np, F);
                status |= _gr_vec_scatter(W.Fim, W.ppos, W.Pim, np, F);
            }

            /* the composites f(d) f(n/d), d the smallest prime factor:
               vectorized when every cofactor n/d <= hi/3 precedes the
               block, otherwise one at a time in increasing order */
            if (hi / 3 < lo)
            {
                for (i = 0; i < nc; i++)
                {
                    W.ia[i] = spf[W.cpos[i]] - 1;
                    W.ib[i] = cof[W.cpos[i]] - 1;
                }
                status |= _gr_vec_gather(W.Ar, tab, W.ia, nc, F);
                status |= _gr_vec_gather(W.Ai, ENTRY(tab, 1, F), W.ia, nc, F);
                status |= _gr_vec_gather(W.Br, tab, W.ib, nc, F);
                status |= _gr_vec_gather(W.Bi, ENTRY(tab, 1, F), W.ib, nc, F);
                status |= _gr_vec_mul(W.T1, W.Ar, W.Br, nc, F);
                status |= _gr_vec_mul(W.T2, W.Ai, W.Bi, nc, F);
                status |= _gr_vec_mul(W.T3, W.Ar, W.Bi, nc, F);
                status |= _gr_vec_mul(W.T4, W.Ai, W.Br, nc, F);
                status |= _gr_vec_sub(W.T1, W.T1, W.T2, nc, F);
                status |= _gr_vec_add(W.T3, W.T3, W.T4, nc, F);
                status |= _gr_vec_scatter(W.Fre, W.cpos, W.T1, nc, F);
                status |= _gr_vec_scatter(W.Fim, W.cpos, W.T3, nc, F);
            }
            else
            {
                /* the primes of the block go to the table first */
                for (i = 0; i < np; i++)
                {
                    k = W.pl[i];
                    if (k <= Mt)
                    {
                        status |= gr_set(ENTRY(tab, k - 1, F), ENTRY(W.Fre, W.ppos[i], F), F);
                        status |= gr_set(ENTRY(tab, k, F), ENTRY(W.Fim, W.ppos[i], F), F);
                    }
                }
                if (lo == 1 && Mt >= 1)
                {
                    status |= gr_one(tab, F);
                    status |= gr_zero(ENTRY(tab, 1, F), F);
                }
                for (i = 0; i < nc; i++)
                {
                    ulong n = lo + 2 * W.cpos[i], d = spf[W.cpos[i]], e = cof[W.cpos[i]];
                    status |= powsum_cmul(ENTRY(W.Fre, W.cpos[i], F), ENTRY(W.Fim, W.cpos[i], F),
                        ENTRY(tab, d - 1, F), ENTRY(tab, d, F), ENTRY(tab, e - 1, F), ENTRY(tab, e, F),
                        s1, s2, F);
                    if (n <= Mt)
                    {
                        status |= gr_set(ENTRY(tab, n - 1, F), ENTRY(W.Fre, W.cpos[i], F), F);
                        status |= gr_set(ENTRY(tab, n, F), ENTRY(W.Fim, W.cpos[i], F), F);
                    }
                }
            }

            /* the table */
            if (lo <= Mt)
            {
                slong m = FLINT_MIN(cnt, (slong) ((Mt - lo) / 2 + 1));
                status |= _gr_vec_scatter(ENTRY(tab, lo - 1, F), W.even, W.Fre, m, F);
                status |= _gr_vec_scatter(ENTRY(tab, lo, F), W.even, W.Fim, m, F);
            }

            /* the partial sums over the pieces of constant floor(log2(N/k)) */
            for (i = 0; i < cnt; )
            {
                slong jj, last;
                ulong bound;
                k = lo + 2 * i;
                jj = FLINT_BIT_COUNT(N / k) - 1;     /* floor(log2(N/k)) */
                bound = N >> jj;                     /* k' <= bound iff J(k') >= jj */
                last = FLINT_MIN(cnt - 1, (slong) ((bound - lo) / 2));
                status |= powsum_tree_sum(s1, ENTRY(W.Fre, i, F), last - i + 1, F);
                status |= powsum_tree_sum(s2, ENTRY(W.Fim, i, F), last - i + 1, F);
                status |= gr_set_other(t1, s1, F, T);
                status |= gr_set_other(t2, s2, F, T);
                status |= gr_add(ENTRY(A->ure, jj, T), ENTRY(A->ure, jj, T), t1, T);
                status |= gr_add(ENTRY(A->uim, jj, T), ENTRY(A->uim, jj, T), t2, T);
                i = last + 1;
            }
        }
    }

    GR_TMP_CLEAR2(t1, t2, T);
    GR_TMP_CLEAR2(s1, s2, F);
    powsum_work_clear(&W);
    *eps_out = eps;
    return status;
}

/* the chunk i of the second phase */
static void
powsum_chunk_worker(slong i, void * args)
{
    const powsum_args_struct * A = args;
    ulong lo0 = A->c0 + 2 * (ulong) i * A->clen;
    ulong hi0 = FLINT_MIN(lo0 + 2 * (A->clen - 1), A->cN);
    A->status[i] = gr_zero(ENTRY(A->Rre, i, A->T), A->T);
    A->status[i] |= gr_zero(ENTRY(A->Rim, i, A->T), A->T);
    A->status[i] |= powsum_range(A, lo0, hi0, ENTRY(A->Rre, i, A->T), ENTRY(A->Rim, i, A->T), A->eps + i);
}

int
gr_powsum_sieved(acb_t res, const acb_t s, ulong N, gr_ctx_t T, gr_ctx_t P,
    gr_ctx_t F_plain, double dm, double da, ulong M)
{
    gr_ctx_struct * F;                  /* the ring of the composites and block sums */
    gr_ptr tab;                         /* f(k) for odd k <= Mt, re at k - 1, im at k */
    gr_ptr ure, uim;                    /* the bucket sums (T) */
    gr_ptr Gre = NULL, Gim = NULL;      /* the table G (T) */
    gr_ptr Gmre = NULL, Gmim = NULL;    /* its midpoints (F; plain) */
    gr_ptr Rre, Rim;                    /* the sum over the roots (T) */
    gr_ptr negsigma, twopi, x2re, x2im, zr, zi, ar, ai, t1, t2, t3, t4;
    gr_ptr tau;
    powsum_work_struct W;
    powsum_args_struct A;
    ulong * small;
    slong * row = NULL;
    slong nsmall, tablen, blk, seg, J, j, i, Glen = 0, nchunks = 0;
    ulong k, sq, Mt, L1, Z = 0;
    arb_t v, sigma_lo;
    gr_ctx_t actx;
    slong wp;
    double eps = 0.0, rG = 0.0;
    int plain = (F_plain != NULL), critical;
    int status = GR_SUCCESS;

    if (N <= 1)
    {
        acb_set_ui(res, N);
        return GR_SUCCESS;
    }
#if FLINT_BITS == 64
    if (N >= (UWORD(1) << 50))
        return GR_UNABLE;
#endif

    F = plain ? F_plain : T;

    /* the table bound: all of it (N/3) for M >= N/3 (and by default,
       M = 0, for N < 2048, where the second phase does not pay),
       otherwise (and by default) at least N^(2/3), so that the table
       G, of about Z^2 / (2 log Z) entries with Z = N/M, is no larger
       than about M */
    if (M >= N / 3 || (M == 0 && N < 2048))
    {
        Mt = N / 3;
        L1 = N;
    }
    else
    {
        M = FLINT_MAX(M, powsum_default_table_bound(N));
        if (M >= N / 3)
        {
            Mt = N / 3;
            L1 = N;
        }
        else
        {
            Mt = M;
            L1 = M;
            Z = N / (M + 1);
        }
    }

    GR_TMP_INIT5(negsigma, twopi, x2re, x2im, zr, T);
    GR_TMP_INIT5(zi, ar, ai, t1, t2, T);
    GR_TMP_INIT4(t3, t4, Rre, Rim, T);
    GR_TMP_INIT(tau, P);

    /* the constants */
    wp = 400;
    gr_ctx_init_real_arb(actx, wp);
    arb_init(v);
    arb_init(sigma_lo);
    arb_const_pi(v, wp);
    arb_mul_2exp_si(v, v, 1);
    status |= gr_set_other(twopi, v, actx, T);
    arb_div(v, acb_imagref(s), v, wp);
    status |= gr_set_other(tau, v, actx, P);
    arb_neg(v, acb_realref(s));
    status |= gr_set_other(negsigma, v, actx, T);
    status |= gr_set_d(t1, -0.5, T);
    critical = (status == GR_SUCCESS && gr_equal(negsigma, t1, T) == T_TRUE);
    {
        arf_t lb;
        arf_init(lb);
        arb_get_lbound_arf(lb, acb_realref(s), 64);
        arb_set_arf(sigma_lo, lb);
        arf_clear(lb);
    }
    arb_clear(v);
    if (status != GR_SUCCESS)
        goto cleanup_early;

    /* the odd primes up to sqrt(N) */
    sq = n_sqrt(N) + 1;
    small = flint_malloc(sizeof(ulong) * (sq / 2 + 2));
    nsmall = 0;
    {
        char * comp = flint_calloc(sq + 1, 1);
        for (k = 3; k <= sq; k += 2)
            if (!comp[k])
            {
                small[nsmall++] = k;
                for (i = k * k; (ulong) i <= sq; i += 2 * k)
                    comp[i] = 1;
            }
        flint_free(comp);
    }

    tablen = (Mt + 1) / 2 + 1;
    tab = gr_heap_init_vec(2 * tablen, F);
    /* the block length: no longer than the number of odd k <= N, so
       that small N do not pay for the full-size scratch; the sieve
       segment: a multiple of it, about sqrt(N)/2 odd numbers (up to
       128 blocks) */
    blk = FLINT_MIN(POWSUM_BLOCK, (slong) ((N + 1) / 2));
    seg = blk * FLINT_MAX(1, FLINT_MIN(128, (slong) (sq / (2 * blk))));

    J = FLINT_BIT_COUNT(N) - 1;         /* floor(log2 N) */
    ure = gr_heap_init_vec(J + 1, T);
    uim = gr_heap_init_vec(J + 1, T);
    status |= gr_zero(Rre, T);
    status |= gr_zero(Rim, T);

    /* f(2), and the table G for the roots beyond L1 */
    powsum_work_init(&W, T, P, F, blk, blk, 0, critical);
    {
        ulong two = 2;
        status |= powsum_prime_terms(x2re, x2im, &two, 1, tau, negsigma, twopi, &W.S);
    }
    if (Z != 0)
    {
        Glen = powsum_G_length(Z);
        Gre = gr_heap_init_vec(Glen, T);
        Gim = gr_heap_init_vec(Glen, T);
        row = flint_calloc(Z + 1, sizeof(slong));
        status |= powsum_G_table(Gre, Gim, row, Z, W.Pre, W.Pim, W.pl, blk, tau, negsigma, twopi, &W.S);
        if (plain)
        {
            /* the midpoints in F, and the largest radius */
            Gmre = gr_heap_init_vec(Glen, F);
            Gmim = gr_heap_init_vec(Glen, F);
            for (i = 0; i < Glen && status == GR_SUCCESS; i++)
            {
                double r1, r2;
                status |= gr_get_interval_mid_rad(t1, t2, ENTRY(Gre, i, T), T);
                status |= gr_set_other(ENTRY(Gmre, i, F), t1, T, F);
                status |= gr_get_d(&r1, t2, T);
                status |= gr_get_interval_mid_rad(t1, t2, ENTRY(Gim, i, T), T);
                status |= gr_set_other(ENTRY(Gmim, i, F), t1, T, F);
                status |= gr_get_d(&r2, t2, T);
                rG = FLINT_MAX(rG, FLINT_MAX(r1, r2));
            }
            rG = rG * (1.0 + 0x1p-40);
            if (!(rG <= 1e-3))
                status = GR_UNABLE;
        }
    }
    powsum_work_clear(&W);

    A.T = T; A.P = P; A.F = F;
    A.plain = plain; A.critical = critical;
    A.N = N; A.Mt = Mt; A.L1 = L1; A.Z = Z;
    A.blk = blk; A.seg = seg; A.nsmall = nsmall; A.small = small;
    A.tab = tab;
    A.Gre = plain ? Gmre : Gre;
    A.Gim = plain ? Gmim : Gim;
    A.row = row;
    A.tau = tau; A.negsigma = negsigma; A.twopi = twopi;
    A.ure = ure; A.uim = uim;
    A.Rre = A.Rim = NULL;
    A.eps = NULL;
    A.status = NULL;

    /* the first phase: the odd k <= L1, in order (the table) */
    if (status == GR_SUCCESS)
        status |= powsum_range(&A, 1, L1 - (L1 % 2 == 0), Rre, Rim, &eps);

    /* the second phase: the odd c in (L1, N], in chunks of a fixed
       length (so that the result does not depend on the number of
       threads), in parallel */
    if (Z != 0 && status == GR_SUCCESS)
    {
        ulong nodd;
        A.c0 = L1 + 1 + (L1 % 2 == 1);
        A.cN = N - (N % 2 == 0);
        nodd = (A.cN - A.c0) / 2 + 1;
        A.clen = FLINT_MAX(UWORD(1) << 18, nodd / 65536 + 1);
        A.clen = FLINT_MAX(A.clen, 16 * (ulong) seg);
        A.clen = ((A.clen + seg - 1) / seg) * seg;
        nchunks = (nodd + A.clen - 1) / A.clen;
        A.Rre = gr_heap_init_vec(nchunks, T);
        A.Rim = gr_heap_init_vec(nchunks, T);
        A.eps = flint_calloc(nchunks, sizeof(double));
        A.status = flint_calloc(nchunks, sizeof(int));
        flint_parallel_do(powsum_chunk_worker, &A, nchunks, 0, FLINT_PARALLEL_STRIDED);
        for (i = 0; i < nchunks; i++)
        {
            status |= A.status[i];
            status |= gr_add(Rre, Rre, ENTRY(A.Rre, i, T), T);
            status |= gr_add(Rim, Rim, ENTRY(A.Rim, i, T), T);
            eps = FLINT_MAX(eps, A.eps[i]);
        }
        gr_heap_clear_vec(A.Rre, nchunks, T);
        gr_heap_clear_vec(A.Rim, nchunks, T);
        flint_free(A.eps);
        flint_free(A.status);
    }

    /* the a priori part of the error */
    if (plain && status == GR_SUCCESS)
    {
        double mu, rho, x, Om, Dt;
        mu = 1.4142135623730951 * (dm + da + dm * da) * (1.0 + 0x1p-40);
        Om = floor(log((double) N) / log(3.0)) + 1.0;     /* >= Omega(n), n odd <= N */
        x = (Om * eps + (Om - 1.0) * mu) * (1.0 + 0x1p-40);
        rho = x * (1.0 + x) * (1.0 + 0x1p-40);          /* >= e^x - 1 for x <= 1 */
        Dt = POWSUM_LOG_BLOCK * da * (1.0 + 2 * POWSUM_LOG_BLOCK * da);
        if (!(eps <= 1e-3) || !(rho <= 1e-3))
            status = GR_UNABLE;
        for (j = 0; j <= J && status == GR_SUCCESS; j++)
        {
            ulong a = (N >> (j + 1)) + 1, b = FLINT_MIN(N >> j, L1);
            double w, e;
            if (a % 2 == 0) a++;
            if (b % 2 == 0) b--;
            w = powsum_weight_bound(a, b, sigma_lo);
            e = (rho + Dt * (1.0 + rho)) * w;
            e = e * (1.0 + 0x1p-40);
            status |= gr_set_d(t1, e, T);
            status |= gr_set_interval_mid_rad(ENTRY(ure, j, T), ENTRY(ure, j, T), t1, T);
            status |= gr_set_interval_mid_rad(ENTRY(uim, j, T), ENTRY(uim, j, T), t1, T);
        }

        /* the roots: the term f(c) G(N/c, d) is computed as F(c) Gm
           (a plain product) with |F(c) - f(c)| <= rho g(c),
           |Gm - G| <= rG and |G| <= Wf = the sum of g(m) over the same
           m as G, so its error is at most
           (rho + mu (1 + rho)) g(c) (Wf + rG) + rG g(c) and its size at
           most (1 + rho) (1 + mu) g(c) (Wf + rG).  As g is completely
           multiplicative, the sum over the roots of g(c) Wf is
           A = sum of g(n) over the n <= N whose odd part is > L1
           (the decomposition n = m c is unique), and the sum of g(c)
           is at most B = the sum of g(k) over the odd k in (L1, N]. */
        if (Z != 0 && status == GR_SUCCESS)
        {
            double A = 0.0, B, X, e;
            ulong a0 = L1 + 1 + (L1 % 2 == 1), b;
            slong aa;

            for (aa = 0; aa <= J && (N >> aa) >= a0; aa++)
            {
                mag_t g2;
                b = N >> aa;
                if (b % 2 == 0) b--;
                if (b < a0)
                    continue;
                mag_init(g2);
                powsum_g_add(g2, UWORD(1) << aa, arb_midref(sigma_lo));
                A += powsum_mag_get_d_up(g2) * powsum_weight_bound(a0, b, sigma_lo) * (1.0 + 0x1p-50);
                mag_clear(g2);
            }
            A *= (1.0 + 0x1p-40);
            B = powsum_weight_bound(a0, N - (N % 2 == 0), sigma_lo);
            X = (A + rG * B) * (1.0 + 0x1p-50);
            e = (rho + mu * (1.0 + rho) + Dt * (1.0 + rho) * (1.0 + mu)) * X + rG * B;
            e = e * (1.0 + 0x1p-40);
            status |= gr_set_d(t1, e, T);
            status |= gr_set_interval_mid_rad(Rre, Rre, t1, T);
            status |= gr_set_interval_mid_rad(Rim, Rim, t1, T);
        }
    }

    /* sum_j f(2)^j U_j with U_j = sum_{m >= j} u_m, by Horner */
    status |= gr_zero(zr, T);
    status |= gr_zero(zi, T);
    status |= gr_zero(ar, T);
    status |= gr_zero(ai, T);
    for (j = J; j >= 0 && status == GR_SUCCESS; j--)
    {
        status |= gr_add(ar, ar, ENTRY(ure, j, T), T);
        status |= gr_add(ai, ai, ENTRY(uim, j, T), T);
        status |= powsum_cmul(t1, t2, zr, zi, x2re, x2im, t3, t4, T);
        status |= gr_add(zr, t1, ar, T);
        status |= gr_add(zi, t2, ai, T);
    }
    status |= gr_add(zr, zr, Rre, T);
    status |= gr_add(zi, zi, Rim, T);
    if (status == GR_SUCCESS)
    {
        status |= gr_set_other(acb_realref(res), zr, T, actx);
        status |= gr_set_other(acb_imagref(res), zi, T, actx);
    }

    flint_free(small);
    gr_heap_clear_vec(tab, 2 * tablen, F);
    if (Z != 0)
    {
        gr_heap_clear_vec(Gre, Glen, T);
        gr_heap_clear_vec(Gim, Glen, T);
        if (plain)
        {
            gr_heap_clear_vec(Gmre, Glen, F);
            gr_heap_clear_vec(Gmim, Glen, F);
        }
        flint_free(row);
    }
    gr_heap_clear_vec(ure, J + 1, T);
    gr_heap_clear_vec(uim, J + 1, T);

cleanup_early:
    arb_clear(sigma_lo);
    gr_ctx_clear(actx);
    GR_TMP_CLEAR5(negsigma, twopi, x2re, x2im, zr, T);
    GR_TMP_CLEAR5(zi, ar, ai, t1, t2, T);
    GR_TMP_CLEAR4(t3, t4, Rre, Rim, T);
    GR_TMP_CLEAR(tau, P);
    return status;
}

/* The dfloat instantiations: terms in dTb balls, phases in dPb balls,
   and (plain) the composites and block sums in plain dT numbers for
   T <= 2, whose weak renormalisation already gives canonical results,
   with the relative error bounds u, u (T = 1) and 8 u^2, 3 u^2 (T = 2)
   of the product and the sum (dev/dfloat_exp_bound.py). */
static int
_dfloat_powsum_sieved_ctx(acb_t res, const acb_t s, ulong N, int tn, int pn, int plain, ulong M)
{
    gr_ctx_t T, P, F;
    double dm, da;
    int status;

    if (tn < 1 || tn > 4 || pn < tn || pn > 4)
        return 0;
    if (plain && tn > 2)
        return 0;
    if (!dfloat_is_supported())
        return 0;

    /* a strong (canonicalising) phase ring: the reduction modulo 1
       cancels the integer part, and the rounding of the canonical
       result to T is a truncation */
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(T, tn, DFLOAT_BALL));
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(P, pn, DFLOAT_BALL | DFLOAT_STRONG));
    if (plain)
    {
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(F, tn, 0));
        if (tn == 1)
        {
            dm = 0x1p-53 * (1.0 + 0x1p-40);
            da = 0x1p-53 * (1.0 + 0x1p-40);
        }
        else
        {
            dm = 8.0 * 0x1p-106 * (1.0 + 0x1p-40);
            da = 3.0 * 0x1p-106 * (1.0 + 0x1p-40);
        }
        status = gr_powsum_sieved(res, s, N, T, P, F, dm, da, M);
        gr_ctx_clear(F);
    }
    else
    {
        status = gr_powsum_sieved(res, s, N, T, P, NULL, 0.0, 0.0, M);
    }
    gr_ctx_clear(T);
    gr_ctx_clear(P);
    return status == GR_SUCCESS;
}

/* the version with all terms in balls (for comparison; the only one
   for T >= 3) */
int
_dfloat_powsum_sieved_ball(acb_t res, const acb_t s, ulong N, int tn, int pn, ulong M)
{
    return _dfloat_powsum_sieved_ctx(res, s, N, tn, pn, 0, M);
}

int
_dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, int tn, int pn, ulong M)
{
    return _dfloat_powsum_sieved_ctx(res, s, N, tn, pn, tn <= 2, M);
}

/* The term precision T (in doubles) so that the result carries about
   prec bits after the binary point: the sum of N terms of size up to
   n^(-sigma) accumulates about 2 sqrt(N) u^T, so the terms need about
   need = prec + log2(N)/2 + 8 bits; the phase precision P so that
   t log N (up to 2^L) keeps need + 6 bits after the binary point (not
   necessarily all of T: at the height of the 10^15-th zero, 132 bits
   need (T, P) = (3, 4), where T full doubles would ask for P = 5). */
int
dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, slong prec)
{
    double t, L;
    int tn, pn;
    slong need;
    ulong M;

    if (N <= 1)
    {
        acb_set_ui(res, N);
        return 1;
    }
    if (!acb_is_finite(s) || arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 10) > 0)
        return 0;
    t = fabs(arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_UP));
    if (!(t < 1e60))
        return 0;
    need = prec + FLINT_BIT_COUNT(N) / 2 + 8;
    tn = (need + 52) / 53;
    tn = FLINT_MAX(tn, 1);
    L = log2(t / (2 * 3.141592653589793) * log((double) N) + 1.0);
    pn = (int) ceil((need + L + 6) / 53);
    /* at the top (P = 4), accept a phase short of the ideal by up to
       8 bits (the result then has up to about 6 bits less than asked,
       as at the height of the 10^17-th zero with 135 bits asked) */
    if (pn == 5 && need + L - 2 <= 4 * 53)
        pn = 4;
    pn = FLINT_MAX(pn, tn);
    if (tn > 3 || pn > 4)
        return 0;
    /* memory: the table of f(k) for odd k <= M (complex, plain for
       T <= 2 and balls for T = 3), with M = N^(2/3) (the default of
       gr_powsum_sieved), and the table G of about Z^2 / (2 log Z)
       complex balls and (for T <= 2) their plain midpoints, Z = N/M */
    {
        double bytes = 2 * 8 * (tn + (tn == 3)), mem, Mt, Z;
        Mt = (N < 2048) ? (double) N / 3 : FLINT_MIN((double) powsum_default_table_bound(N), (double) N / 3);
        mem = (Mt / 2 + 1) * bytes;
        if (Mt < (double) N / 3)
        {
            Z = (double) N / Mt;
            mem += (Z + 1 + Z * Z / (2 * log(Z))) * (2 * 8 * (tn + 1) + (tn <= 2 ? bytes : 0));
        }
        if (mem > DFLOAT_POWSUM_MAX_BYTES)
            return 0;
        M = 0;
    }
    return _dfloat_powsum_sieved(res, s, N, tn, pn, M);
}

/* ------------------------------------------------------------------------
   The sums over j of Platt's multi-evaluation (acb_dirichlet_platt_multieval,
   _platt_smk): for the buckets m of the j with smk_points[m] <= j <
   smk_points[m + 1] (1 <= j <= J), S[k][m] = sum_j z_j b_j^k (0 <= k < K),
   where z_j = j^(-1/2) exp(-i t0 log(j sqrt(pi))) = c f(j),
   f(j) = j^(-s) with s = 1/2 + i t0, c = exp(-i t0 log sqrt(pi)), and
   b_j = log(j sqrt(pi)) / (2 pi) - m/B, |b_j| <= 1/(2B).

   The terms: f is completely multiplicative and log additive, so for a
   composite j with smallest prime factor d and cofactor j/d <= M,
   f(j) = f(d) f(j/d) and log j = log d + log(j/d) from a table (d3b
   balls, k <= M); 1, the primes and the other composites are computed
   as the prime terms of the power sum (phases in dPb balls, the rest in
   d3b balls); b_j in d3b balls.  The j are processed in dyadic ranges
   while 2^i <= M (the cofactors of a range lie in the previous ones),
   each in chunks of a fixed length, in parallel.

   The rest of the multi-evaluation amplifies an error in S[k] by about
   2^(24 + 8.6 k) (heuristic parameters, B = 4096) while |S[k]| is about
   W 2^(-13 k) (W the sum of j^(-1/2) over the bucket), so for a grid
   good to about 75 bits after the binary point S[0] needs
   triple-double accuracy, S[1] and S[2] double-double with sums of
   logarithmic depth, beyond k ~ 13 double: S[0] is summed in d3b
   balls; the moments 1 <= k < PSMK_KD from p_j = z_j (the midpoints
   rounded to double-words) in plain double-double arithmetic and
   beyond in plain double arithmetic, one k at a time over all j of a
   block, with a priori error bounds: with Z_j = |Re p~| + |Im p~| and
   R_j the radii of z (both parts, including the rounding to
   double-words), rho and beta bounds for the radius and the magnitude
   of b (and of its midpoint), the computed sum over a bucket of the
   moment k differs from the exact one by at most, in each part,

     beta^k sum R_j + sum Z_j (k rho beta^(k-1) + beta^k (e_k + s_k (1 + e_k)))

   where e_k bounds the relative error of the computed product (k
   double-double multiplications DWTimesDW1 of relative error at most
   7u^2, then a rounding to double and double multiplications of
   relative error u) and s_k the summation error (D times 4u^2 for
   AccurateDWPlusDW, of relative error at most 3u^2/(1 - 4u), or D u
   in double, D the depth of the summation tree of the bucket), plus
   an allowance of 2^-1060 per operation for underflow. */

#define PSMK_KB 1       /* the moment k = 0 in d3b balls */
#define PSMK_KD 14      /* then k < PSMK_KD in double-double, then double */
/* the table of f(k), log k, k <= M */
#if FLINT_BITS == 64
#define PSMK_TABLE_BYTES 5e8
#else
#define PSMK_TABLE_BYTES 1e8
#endif

typedef struct
{
    double hi, lo;
}
psmk_dd;

/* AccurateDWPlusDW (Joldes, Muller, Popescu 2017): relative error at
   most 3u^2/(1 - 4u) of |a + b| */
FLINT_FORCE_INLINE psmk_dd psmk_dd_add(psmk_dd a, psmk_dd b)
{
    psmk_dd r;
    double s, e, t, f;
    DFLOAT_TWO_SUM(s, e, a.hi, b.hi);
    DFLOAT_TWO_SUM(t, f, a.lo, b.lo);
    e += t;
    DFLOAT_FAST_TWO_SUM(s, e, s, e);
    e += f;
    DFLOAT_FAST_TWO_SUM(r.hi, r.lo, s, e);
    return r;
}

static void
psmk_arb_set_dd(arb_t x, double hi, double lo)
{
    arf_t t;
    arf_init(t);
    arf_set_d(arb_midref(x), hi);
    arf_set_d(t, lo);
    arf_add(arb_midref(x), arb_midref(x), t, ARF_PREC_EXACT, ARF_RND_DOWN);
    mag_zero(arb_radref(x));
    arf_clear(t);
}

/* Compact complex balls: re hi, re lo, im hi, im lo (double-words, or
   any doubles), and a radius for both parts. */
static void
psmk_acb_from5(acb_t x, const double * e)
{
    mag_t r;
    psmk_arb_set_dd(acb_realref(x), e[0], e[1]);
    psmk_arb_set_dd(acb_imagref(x), e[2], e[3]);
    mag_init(r);
    mag_set_d(r, e[4]);
    mag_set(arb_radref(acb_realref(x)), r);
    mag_set(arb_radref(acb_imagref(x)), r);
    mag_clear(r);
}

/* the midpoint rounded to a double-word, the rounding error added to
   the radius */
static double
psmk_arb_to_dd(double * hi, double * lo, const arb_t x)
{
    arf_t t;
    mag_t m;
    double r;
    arf_init(t);
    mag_init(m);
    *hi = arf_get_d(arb_midref(x), ARF_RND_NEAR);
    arf_set_d(t, *hi);
    arf_sub(t, arb_midref(x), t, ARF_PREC_EXACT, ARF_RND_DOWN);
    *lo = arf_get_d(t, ARF_RND_NEAR);
    {
        arf_t u;
        arf_init(u);
        arf_set_d(u, *lo);
        arf_sub(t, t, u, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_clear(u);
    }
    arf_get_mag(m, t);
    mag_add(m, m, arb_radref(x));
    r = mag_get_d(m);
    arf_clear(t);
    mag_clear(m);
    return r;
}

static void
psmk_acb_to5(double * e, const acb_t x)
{
    double r1, r2;
    r1 = psmk_arb_to_dd(e, e + 1, acb_realref(x));
    r2 = psmk_arb_to_dd(e + 2, e + 3, acb_imagref(x));
    e[4] = FLINT_MAX(r1, r2);
}

/* x *= b elementwise (double-double, DWTimesDW1), for the real and
   imaginary parts; written for vectorization */
static void
psmk_dd_vec_mul(double * restrict xh, double * restrict xl,
    double * restrict yh, double * restrict yl,
    const double * restrict bh, const double * restrict bl, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        double p, e, t;
        DFLOAT_TWO_PROD(p, e, xh[i], bh[i]);
        t = xh[i] * bl[i] + xl[i] * bh[i];
        e = e + t;
        DFLOAT_FAST_TWO_SUM(xh[i], xl[i], p, e);
        DFLOAT_TWO_PROD(p, e, yh[i], bh[i]);
        t = yh[i] * bl[i] + yl[i] * bh[i];
        e = e + t;
        DFLOAT_FAST_TWO_SUM(yh[i], yl[i], p, e);
    }
}

#define PSMK_LANES 8
#define PSMK_RUN 64

/* The sum of x[0..len) (double-double, AccurateDWPlusDW): runs of
   PSMK_RUN elements summed with PSMK_LANES interleaved accumulators
   (to break the dependency chains) and combined pairwise, then the run
   sums pairwise.  Every element takes part in at most *depth
   additions, so the error is at most depth 4u^2 (1 + 4u^2)^depth times
   the sum of |x|. */
static psmk_dd
psmk_dd_vec_sum(const double * restrict xh, const double * restrict xl,
    slong len, slong * depth)
{
    psmk_dd runs[4096 / PSMK_RUN + 2];
    slong nr = 0, i0, i, l, n, d;

    for (i0 = 0; i0 < len; i0 += PSMK_RUN)
    {
        double ah[PSMK_LANES], al[PSMK_LANES];
        slong i1 = FLINT_MIN(len, i0 + PSMK_RUN);

        for (l = 0; l < PSMK_LANES; l++)
            ah[l] = al[l] = 0.0;
        for (i = i0; i + PSMK_LANES <= i1; i += PSMK_LANES)
        {
            for (l = 0; l < PSMK_LANES; l++)
            {
                double s, e, t, f;
                DFLOAT_TWO_SUM(s, e, ah[l], xh[i + l]);
                DFLOAT_TWO_SUM(t, f, al[l], xl[i + l]);
                e += t;
                DFLOAT_FAST_TWO_SUM(s, e, s, e);
                e += f;
                DFLOAT_FAST_TWO_SUM(ah[l], al[l], s, e);
            }
        }
        for (l = 0; i < i1; i++, l++)
        {
            double s, e, t, f;
            DFLOAT_TWO_SUM(s, e, ah[l], xh[i]);
            DFLOAT_TWO_SUM(t, f, al[l], xl[i]);
            e += t;
            DFLOAT_FAST_TWO_SUM(s, e, s, e);
            e += f;
            DFLOAT_FAST_TWO_SUM(ah[l], al[l], s, e);
        }
        for (n = PSMK_LANES / 2; n >= 1; n /= 2)
        {
            for (l = 0; l < n; l++)
            {
                psmk_dd x, y;
                x.hi = ah[l]; x.lo = al[l];
                y.hi = ah[l + n]; y.lo = al[l + n];
                x = psmk_dd_add(x, y);
                ah[l] = x.hi; al[l] = x.lo;
            }
        }
        runs[nr].hi = ah[0];
        runs[nr].lo = al[0];
        nr++;
    }

    /* depth within a run: PSMK_RUN / PSMK_LANES + log2(PSMK_LANES) */
    d = PSMK_RUN / PSMK_LANES + 3;
    while (nr > 1)
    {
        for (i = 0; i < nr / 2; i++)
            runs[i] = psmk_dd_add(runs[2 * i], runs[2 * i + 1]);
        if (nr % 2)
            runs[nr / 2] = runs[nr - 1];
        nr = (nr + 1) / 2;
        d++;
    }
    *depth = d;
    if (len == 0)
    {
        runs[0].hi = runs[0].lo = 0.0;
    }
    return runs[0];
}

/* the same in double arithmetic (error at most depth u (1 + u)^depth
   times the sum of |x|) */
static double
psmk_d_vec_sum(const double * restrict x, slong len, slong * depth)
{
    double runs[4096 / PSMK_RUN + 2];
    slong nr = 0, i0, i, l, n, d;

    for (i0 = 0; i0 < len; i0 += PSMK_RUN)
    {
        double a[PSMK_LANES];
        slong i1 = FLINT_MIN(len, i0 + PSMK_RUN);
        for (l = 0; l < PSMK_LANES; l++)
            a[l] = 0.0;
        for (i = i0; i + PSMK_LANES <= i1; i += PSMK_LANES)
            for (l = 0; l < PSMK_LANES; l++)
                a[l] += x[i + l];
        for (l = 0; i < i1; i++, l++)
            a[l] += x[i];
        for (n = PSMK_LANES / 2; n >= 1; n /= 2)
            for (l = 0; l < n; l++)
                a[l] += a[l + n];
        runs[nr++] = a[0];
    }
    d = PSMK_RUN / PSMK_LANES + 3;
    while (nr > 1)
    {
        for (i = 0; i < nr / 2; i++)
            runs[i] = runs[2 * i] + runs[2 * i + 1];
        if (nr % 2)
            runs[nr / 2] = runs[nr - 1];
        nr = (nr + 1) / 2;
        d++;
    }
    *depth = d;
    return (len == 0) ? 0.0 : runs[0];
}

typedef struct
{
    /* shared */
    const ulong * pts;          /* bucket thresholds, pts[N] = UWORD_MAX */
    slong A, B, K, N, pn, prec;
    ulong J;
    arb_srcptr t0;
    /* the chunks: [cstart[c], cend[c]]; this call does the chunks
       c0 + c of the current range */
    const ulong * cstart, * cend;
    slong c0;
    /* the sieve: the odd primes up to sqrt(J), and the table of f(k)
       and log k for k <= M (d3b balls; written by the chunks of the
       ranges below M, read by the later ones) */
    const ulong * small;
    slong nsmall;
    ulong M;
    gr_ptr tre, tim, tlg;
    /* per chunk output: the buckets m0..m1 and their sums, K per bucket */
    slong * m0, * m1;
    acb_ptr * res;
    int * status;
}
psmk_args;

static void
psmk_chunk(slong cc, void * arg_ptr)
{
    const psmk_args * a = arg_ptr;
    slong c = a->c0 + cc;
    gr_ctx_t T, P;
    powsum_scratch_struct S;
    gr_ptr re, im, bas, lt, mv, SbR, SbI, s1;
    gr_ptr negsigma, twopi, tau, inv2pi, c0;
    ulong * list, jlo, jhi, jstart, jend;
    slong * mj;
    slong blk = POWSUM_BLOCK, K = a->K, B = a->B, N = a->N, i, k, m, mfirst, nb, q;
    slong KB = FLINT_MIN(PSMK_KB, K);
    double * sums, * sumZ, * sumR, rho = 0.0, beta = 0.0, * Lb, * nadd;
    double * prh, * prl, * pih, * pil, * bh, * bl, * dre, * dim, * db;
    ulong * spf, * snext, * dlist, * fnext = NULL, * fr = NULL, * fx = NULL, * fy = NULL;
    int split;
    slong * dpos, * cpos, * ia, * ib, nsmall = a->nsmall, jtop = 0;
    gr_ptr Ar, Ai, Al, Br, Bi, Bl, X1, X2;
    acb_ptr res;
    arb_t v, w;
    gr_ctx_t actx;
    int status = GR_SUCCESS;

    jstart = a->cstart[c];
    jend = a->cend[c];
    split = (jstart > a->M);

    GR_MUST_SUCCEED(gr_ctx_init_dfloat(T, 3, DFLOAT_BALL));
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(P, a->pn, DFLOAT_BALL | DFLOAT_STRONG));
    GR_TMP_INIT5(negsigma, twopi, inv2pi, c0, s1, T);
    GR_TMP_INIT(tau, P);
    gr_ctx_init_real_arb(actx, 400);
    arb_init(v);
    arb_init(w);

    arb_const_pi(v, 400);
    arb_mul_2exp_si(v, v, 1);                       /* 2 pi */
    status |= gr_set_other(twopi, v, actx, T);
    arb_div(w, a->t0, v, 400);
    status |= gr_set_other(tau, w, actx, P);        /* t0 / (2 pi) */
    arb_inv(w, v, 400);
    status |= gr_set_other(inv2pi, w, actx, T);
    arb_const_sqrt_pi(w, 400);
    arb_log(w, w, 400);
    arb_div(w, w, v, 400);
    status |= gr_set_other(c0, w, actx, T);         /* log(sqrt(pi)) / (2 pi) */
    status |= gr_set_d(negsigma, -0.5, T);

    S.T = T;
    S.P = P;
    S.critical = 1;
    S.xp = gr_heap_init_vec(blk, P);
    S.lp = gr_heap_init_vec(blk, P);
    S.ph = gr_heap_init_vec(blk, P);
    S.nv = gr_heap_init_vec(blk, P);
    S.fr = gr_heap_init_vec(blk, T);
    S.lt = gr_heap_init_vec(blk, T);
    S.c = gr_heap_init_vec(blk, T);
    S.sn = gr_heap_init_vec(blk, T);
    S.m = gr_heap_init_vec(blk, T);
    re = gr_heap_init_vec(blk, T);
    im = gr_heap_init_vec(blk, T);
    bas = gr_heap_init_vec(blk, T);
    lt = gr_heap_init_vec(blk, T);
    mv = gr_heap_init_vec(blk, T);
    list = flint_malloc(sizeof(ulong) * blk);
    mj = flint_malloc(sizeof(slong) * blk);
    spf = flint_malloc(sizeof(ulong) * blk);
    dlist = flint_malloc(sizeof(ulong) * blk);
    snext = flint_malloc(sizeof(ulong) * (nsmall + 1));
    dpos = flint_malloc(sizeof(slong) * blk);
    cpos = flint_malloc(sizeof(slong) * blk);
    ia = flint_malloc(sizeof(slong) * blk);
    ib = flint_malloc(sizeof(slong) * blk);
    Ar = gr_heap_init_vec(blk, T);
    Ai = gr_heap_init_vec(blk, T);
    Al = gr_heap_init_vec(blk, T);
    Br = gr_heap_init_vec(blk, T);
    Bi = gr_heap_init_vec(blk, T);
    Bl = gr_heap_init_vec(blk, T);
    X1 = gr_heap_init_vec(blk, T);
    X2 = gr_heap_init_vec(blk, T);

    /* beyond M (the last range), the full factorization over the primes
       up to sqrt(J) (the first multiple >= jstart of 2 and of each odd
       prime), for the split j = x y with x, y <= M */
    if (split)
    {
        fnext = flint_malloc(sizeof(ulong) * (nsmall + 1));
        fr = flint_malloc(sizeof(ulong) * blk);
        fx = flint_malloc(sizeof(ulong) * blk);
        fy = flint_malloc(sizeof(ulong) * blk);
        fnext[nsmall] = ((jstart + 1) / 2) * 2;
        for (i = 0; i < nsmall; i++)
        {
            ulong qq = a->small[i];
            fnext[i] = ((jstart + qq - 1) / qq) * qq;
        }
    }

    /* the sieve state: the first multiple >= max(q^2, jstart) of each
       odd prime q (the even j have the smallest prime factor 2) */
    for (i = 0; i < nsmall; i++)
    {
        ulong qq = a->small[i], mm;
        if (qq * qq >= jstart)
            mm = qq * qq;
        else
            mm = ((jstart + qq - 1) / qq) * qq;
        snext[i] = mm;
    }

    /* the buckets of the chunk: the last m with pts[m] <= jstart, and
       with pts[m] <= jend (binary searches) */
    {
        slong lo = 0, hi = N - 1;
        while (lo < hi)
        {
            slong mid = (lo + hi + 1) / 2;
            if (a->pts[mid] <= jstart) lo = mid; else hi = mid - 1;
        }
        mfirst = m = lo;
        hi = N - 1;
        while (lo < hi)
        {
            slong mid = (lo + hi + 1) / 2;
            if (a->pts[mid] <= jend) lo = mid; else hi = mid - 1;
        }
        nb = lo - mfirst + 1;
    }
    sums = flint_calloc(nb * K * 4, sizeof(double));
    sumZ = flint_calloc(nb, sizeof(double));
    sumR = flint_calloc(nb, sizeof(double));
    Lb = flint_calloc(nb, sizeof(double));
    nadd = flint_calloc(nb * K, sizeof(double));
    prh = flint_malloc(sizeof(double) * 9 * blk);
    prl = prh + blk;
    pih = prl + blk;
    pil = pih + blk;
    bh = pil + blk;
    bl = bh + blk;
    dre = bl + blk;
    dim = dre + blk;
    db = dim + blk;
    SbR = gr_heap_init_vec(nb * KB, T);
    SbI = gr_heap_init_vec(nb * KB, T);

    for (jlo = jstart; jlo <= jend && status == GR_SUCCESS; jlo = jhi + 1)
    {
        slong len, seg0;
        const d3b_struct * zr, * zi, * bb;

        jhi = FLINT_MIN(jend, jlo + blk - 1);
        len = jhi - jlo + 1;

        for (i = 0; i < len; i++)
        {
            ulong j = jlo + i;
            list[i] = j;
            while (m < N - 1 && a->pts[m + 1] <= j)
                m++;
            mj[i] = m;
        }

        /* the smallest prime factors (0: j = 1 or a prime) */
        for (i = 0; i < len; i++)
            spf[i] = ((jlo + i) % 2 == 0 && jlo + i > 2) ? 2 : 0;
        while (jtop < nsmall && a->small[jtop] * a->small[jtop] <= jhi)
            jtop++;
        for (k = jtop - 1; k >= 0; k--)
        {
            ulong qq = a->small[k], mm = snext[k];
            for ( ; mm <= jhi; mm += qq)
                if (spf[mm - jlo] != 2)
                    spf[mm - jlo] = qq;
            snext[k] = mm;
        }

        /* beyond M: the factorization j = (prime powers p^e, p <= sqrt(J))
           R, R = 1 or a prime > sqrt(J); the prime powers, in decreasing
           order of p, go to x while x p^e <= M, otherwise to y */
        if (split)
        {
            slong qi;
            for (i = 0; i < len; i++)
            {
                fr[i] = jlo + i;
                fx[i] = fy[i] = 1;
            }
            for (qi = nsmall; qi >= 0; qi--)
            {
                ulong qq = (qi == nsmall) ? 2 : a->small[qi], mm;
                for (mm = fnext[qi]; mm <= jhi; mm += qq)
                {
                    ulong ii = mm - jlo, pe = qq;
                    fr[ii] /= qq;
                    while (fr[ii] % qq == 0)
                    {
                        fr[ii] /= qq;
                        pe *= qq;
                    }
                    if (fx[ii] <= a->M / pe)
                        fx[ii] *= pe;
                    else
                        fy[ii] *= pe;
                }
                fnext[qi] = mm;
            }
        }

        /* f(j) and log(j): for a composite j with smallest prime factor
           d and cofactor j/d <= M, as f(d) f(j/d) and log(d) + log(j/d)
           from the table; beyond M, otherwise, as f(x') f(y') for a
           split j = x' y' with 2 <= x', y' <= M from the factorization
           ((x, y R) or (x R, y)); directly for 1, the primes and the
           rest */
        {
            slong nd = 0, nc = 0, nt = 0;
            for (i = 0; i < len; i++)
            {
                ulong j = jlo + i, d = spf[i], x = 0, y = 0;
                if (d != 0 && j / d <= a->M)
                {
                    x = d;
                    y = j / d;
                }
                else if (split && d != 0)
                {
                    ulong R = fr[i];
                    if (fy[i] <= a->M / R && fx[i] >= 2 && fy[i] * R >= 2)
                    {
                        x = fx[i];
                        y = fy[i] * R;
                    }
                    else if (fx[i] <= a->M / R && fy[i] >= 2 && fx[i] * R >= 2)
                    {
                        x = fx[i] * R;
                        y = fy[i];
                    }
                    /* the identity x y = j is what makes this exact */
                    if (x != 0 && (x > a->M || y > a->M || x * y != j))
                        x = y = 0;
                }

                if (x == 0)
                {
                    dlist[nd] = j;
                    dpos[nd] = i;
                    nd++;
                }
                else
                {
                    ia[nc] = x;
                    ib[nc] = y;
                    cpos[nc] = i;
                    nc++;
                }
                nt += (j <= a->M);
            }

            status |= powsum_prime_terms(X1, X2, dlist, nd, tau, negsigma, twopi, &S);
            status |= _gr_vec_scatter(re, dpos, X1, nd, T);
            status |= _gr_vec_scatter(im, dpos, X2, nd, T);
            status |= _gr_vec_set_other(X1, S.lp, P, nd, T);
            status |= _gr_vec_scatter(lt, dpos, X1, nd, T);

            if (nc > 0)
            {
                status |= _gr_vec_gather(Ar, a->tre, ia, nc, T);
                status |= _gr_vec_gather(Ai, a->tim, ia, nc, T);
                status |= _gr_vec_gather(Al, a->tlg, ia, nc, T);
                status |= _gr_vec_gather(Br, a->tre, ib, nc, T);
                status |= _gr_vec_gather(Bi, a->tim, ib, nc, T);
                status |= _gr_vec_gather(Bl, a->tlg, ib, nc, T);
                status |= _gr_vec_add(Al, Al, Bl, nc, T);
                status |= _gr_vec_scatter(lt, cpos, Al, nc, T);
                status |= _gr_vec_mul(X1, Ar, Br, nc, T);
                status |= _gr_vec_mul(X2, Ai, Bi, nc, T);
                status |= _gr_vec_sub(X1, X1, X2, nc, T);
                status |= _gr_vec_scatter(re, cpos, X1, nc, T);
                status |= _gr_vec_mul(X1, Ar, Bi, nc, T);
                status |= _gr_vec_mul(X2, Ai, Br, nc, T);
                status |= _gr_vec_add(X1, X1, X2, nc, T);
                status |= _gr_vec_scatter(im, cpos, X1, nc, T);
            }

            /* the table (these j, the first nt of the block, are not
               read by this range) */
            if (nt > 0)
            {
                status |= _gr_vec_set(ENTRY(a->tre, jlo, T), re, nt, T);
                status |= _gr_vec_set(ENTRY(a->tim, jlo, T), im, nt, T);
                status |= _gr_vec_set(ENTRY(a->tlg, jlo, T), lt, nt, T);
            }
        }

        /* b_j = log(j) / (2 pi) + c0 - m/B (m/B as a ball, per bucket) */
        status |= _gr_vec_mul_scalar(bas, lt, len, inv2pi, T);
        status |= _gr_vec_add_scalar(bas, bas, len, c0, T);
        {
            slong seg0, seg1;
            for (seg0 = 0; seg0 < len; seg0 = seg1)
            {
                seg1 = seg0 + 1;
                while (seg1 < len && mj[seg1] == mj[seg0])
                    seg1++;
                status |= gr_set_si(mv, mj[seg0], T);
                status |= gr_div_si(mv, mv, B, T);
                status |= _gr_vec_sub_scalar(ENTRY(bas, seg0, T), ENTRY(bas, seg0, T), seg1 - seg0, mv, T);
            }
        }

        /* the moment k = 0, in balls */
        for (seg0 = 0; seg0 < len; )
        {
            slong seg1 = seg0 + 1, qq = mj[seg0] - mfirst;
            while (seg1 < len && mj[seg1] == mj[seg0])
                seg1++;
            status |= _gr_vec_sum(s1, ENTRY(re, seg0, T), seg1 - seg0, T);
            status |= gr_add(ENTRY(SbR, qq, T), ENTRY(SbR, qq, T), s1, T);
            status |= _gr_vec_sum(s1, ENTRY(im, seg0, T), seg1 - seg0, T);
            status |= gr_add(ENTRY(SbI, qq, T), ENTRY(SbI, qq, T), s1, T);
            seg0 = seg1;
        }
        if (K <= 1)
            continue;
        if (status != GR_SUCCESS)
            break;

        /* the moments k >= 1 in plain arithmetic, from p = z: the
           conversions to double-double (TwoSum of the leading
           component and the rounded sum of the other two, an error of
           at most u (|d1| + |d2|)), then, for each k, the next power
           and the bucket sums (psmk_dd_vec_sum), over all j of the
           block */
        zr = (const d3b_struct *) re;
        zi = (const d3b_struct *) im;
        bb = (const d3b_struct *) bas;

        for (i = 0; i < len; i++)
        {
            double Z, R, br;
            slong qq = mj[i] - mfirst;

            DFLOAT_TWO_SUM(prh[i], prl[i], zr[i].d[0], zr[i].d[1] + zr[i].d[2]);
            DFLOAT_TWO_SUM(pih[i], pil[i], zi[i].d[0], zi[i].d[1] + zi[i].d[2]);
            DFLOAT_TWO_SUM(bh[i], bl[i], bb[i].d[0], bb[i].d[1] + bb[i].d[2]);

            Z = fabs(prh[i]) + fabs(prl[i]) + fabs(pih[i]) + fabs(pil[i]);
            R = zr[i].rad + zi[i].rad + (fabs(zr[i].d[1]) + fabs(zr[i].d[2]) + fabs(zi[i].d[1]) + fabs(zi[i].d[2])) * 0x1p-52;
            br = bb[i].rad + (fabs(bb[i].d[1]) + fabs(bb[i].d[2])) * 0x1p-52;
            if (!(R < D_INF) || !(br < D_INF) || !(Z < D_INF))
            {
                status = GR_UNABLE;
                break;
            }
            sumZ[qq] += Z;
            sumR[qq] += R;
            Lb[qq] += 1.0;
            rho = FLINT_MAX(rho, br);
            beta = FLINT_MAX(beta, fabs(bh[i]) + fabs(bl[i]) + br);
        }
        if (status != GR_SUCCESS)
            break;

        for (k = 1; k < K; k++)
        {
            slong seg0, seg1;
            int dd = (k < PSMK_KD);

            /* p = z b^k: k double-double multiplications for k < KD,
               then a rounding to double and double multiplications */
            if (dd)
                psmk_dd_vec_mul(prh, prl, pih, pil, bh, bl, len);
            else
            {
                if (k == PSMK_KD)
                {
                    for (i = 0; i < len; i++)
                    {
                        dre[i] = prh[i] + prl[i];
                        dim[i] = pih[i] + pil[i];
                        db[i] = bh[i] + bl[i];
                    }
                }
                for (i = 0; i < len; i++)
                {
                    dre[i] *= db[i];
                    dim[i] *= db[i];
                }
            }

            for (seg0 = 0; seg0 < len; seg0 = seg1)
            {
                slong qq = mj[seg0] - mfirst, d1, d2;
                double * acc = sums + (qq * K + k) * 4;
                seg1 = seg0 + 1;
                while (seg1 < len && mj[seg1] == mj[seg0])
                    seg1++;

                if (dd)
                {
                    psmk_dd x, y;
                    x = psmk_dd_vec_sum(prh + seg0, prl + seg0, seg1 - seg0, &d1);
                    y.hi = acc[0]; y.lo = acc[1];
                    y = psmk_dd_add(y, x);
                    acc[0] = y.hi; acc[1] = y.lo;
                    x = psmk_dd_vec_sum(pih + seg0, pil + seg0, seg1 - seg0, &d2);
                    y.hi = acc[2]; y.lo = acc[3];
                    y = psmk_dd_add(y, x);
                    acc[2] = y.hi; acc[3] = y.lo;
                }
                else
                {
                    acc[0] += psmk_d_vec_sum(dre + seg0, seg1 - seg0, &d1);
                    acc[2] += psmk_d_vec_sum(dim + seg0, seg1 - seg0, &d2);
                }
                /* the depth of the bucket sum: one more addition */
                nadd[qq * K + k] = FLINT_MAX(nadd[qq * K + k], (double) FLINT_MAX(d1, d2)) + 1.0;
            }
        }
    }

    /* the output: k = 0 from the balls; for k >= 1, the plain sums
       with the error bound (per part)

         beta^k sum R + sum Z (k rho beta^(k-1) + beta^k (e_k + s_k (1 + e_k)))

       where e_k = k 8u^2 for k < KD (k double-double multiplications
       of relative error at most 7u^2), and (KD - 1) 8u^2 + (k - KD + 2) u
       beyond (a rounding to double and k - KD + 1 double
       multiplications), and s_k = D 4u^2 (or D u in double) for a
       bucket sum of depth D (every term takes part in at most D
       additions of relative error at most 3u^2/(1 - 4u), or u); plus
       2^-1060 per operation for underflow. */
    res = _acb_vec_init(nb * K);
    if (status == GR_SUCCESS)
    {
        const double u = 0x1p-53, udd = 8.0 * 0x1p-106, add = 4.0 * 0x1p-106;
        mag_t e;
        mag_init(e);
        rho *= DFLOAT_RAD_SLACK;
        beta *= DFLOAT_RAD_SLACK;
        for (q = 0; q < nb && status == GR_SUCCESS; q++)
        {
            double bk, bk1, sZ = sumZ[q] * DFLOAT_RAD_SLACK, sR = sumR[q] * DFLOAT_RAD_SLACK, L = Lb[q];

            status |= gr_set_other(acb_realref(res + q * K), ENTRY(SbR, q, T), T, actx);
            status |= gr_set_other(acb_imagref(res + q * K), ENTRY(SbI, q, T), T, actx);

            bk1 = 1.0;
            bk = beta;
            for (k = 1; k < K; k++)
            {
                const double * sp = sums + (q * K + k) * 4;
                double ek, sk, err, D = nadd[q * K + k];
                if (k < PSMK_KD)
                {
                    ek = k * udd;
                    sk = D * add * 1.0000001;
                }
                else
                {
                    ek = (PSMK_KD - 1) * udd + (k - PSMK_KD + 2) * u * 1.0000001;
                    sk = D * u * 1.0000001;
                }
                err = bk * sR + sZ * (k * rho * bk1 + bk * (ek + sk * (1.0 + ek)));
                err += L * (k + D + 4) * 0x1p-1060;
                err *= DFLOAT_RAD_SLACK;
                psmk_arb_set_dd(acb_realref(res + q * K + k), sp[0], sp[1]);
                psmk_arb_set_dd(acb_imagref(res + q * K + k), sp[2], sp[3]);
                mag_set_d(e, err);
                arb_add_error_mag(acb_realref(res + q * K + k), e);
                arb_add_error_mag(acb_imagref(res + q * K + k), e);
                bk1 = bk;
                bk = bk * beta * DFLOAT_RAD_SLACK;
            }
        }
        mag_clear(e);
    }

    a->m0[c] = mfirst;
    a->m1[c] = mfirst + nb - 1;
    a->res[c] = res;
    a->status[c] = status;

    flint_free(sums);
    flint_free(sumZ);
    flint_free(sumR);
    flint_free(Lb);
    flint_free(nadd);
    flint_free(prh);
    flint_free(list);
    flint_free(mj);
    flint_free(spf);
    flint_free(dlist);
    flint_free(snext);
    if (split)
    {
        flint_free(fnext);
        flint_free(fr);
        flint_free(fx);
        flint_free(fy);
    }
    flint_free(dpos);
    flint_free(cpos);
    flint_free(ia);
    flint_free(ib);
    gr_heap_clear_vec(Ar, blk, T);
    gr_heap_clear_vec(Ai, blk, T);
    gr_heap_clear_vec(Al, blk, T);
    gr_heap_clear_vec(Br, blk, T);
    gr_heap_clear_vec(Bi, blk, T);
    gr_heap_clear_vec(Bl, blk, T);
    gr_heap_clear_vec(X1, blk, T);
    gr_heap_clear_vec(X2, blk, T);
    gr_heap_clear_vec(SbR, nb * KB, T);
    gr_heap_clear_vec(SbI, nb * KB, T);
    gr_heap_clear_vec(S.xp, blk, P);
    gr_heap_clear_vec(S.lp, blk, P);
    gr_heap_clear_vec(S.ph, blk, P);
    gr_heap_clear_vec(S.nv, blk, P);
    gr_heap_clear_vec(S.fr, blk, T);
    gr_heap_clear_vec(S.lt, blk, T);
    gr_heap_clear_vec(S.c, blk, T);
    gr_heap_clear_vec(S.sn, blk, T);
    gr_heap_clear_vec(S.m, blk, T);
    gr_heap_clear_vec(re, blk, T);
    gr_heap_clear_vec(im, blk, T);
    gr_heap_clear_vec(bas, blk, T);
    gr_heap_clear_vec(lt, blk, T);
    gr_heap_clear_vec(mv, blk, T);
    GR_TMP_CLEAR5(negsigma, twopi, inv2pi, c0, s1, T);
    GR_TMP_CLEAR(tau, P);
    gr_ctx_clear(actx);
    gr_ctx_clear(T);
    gr_ctx_clear(P);
    arb_clear(v);
    arb_clear(w);
}

int
_dfloat_platt_smk_dd(double * S5, const fmpz * smk_points, const arb_t t0,
    slong A, slong B, ulong J, slong K, int pn, ulong M)
{
    psmk_args a;
    slong N = A * B, nchunks, c, m, k, nsmall, nranges, r, * rstart;
    ulong * pts, * small, * cs, * ce;
    gr_ctx_t Tt;
    acb_t y;
    int status = GR_SUCCESS;

    if (pn < 3 || pn > 4 || !dfloat_is_supported() || J < 1)
        return 0;

    pts = flint_malloc(sizeof(ulong) * (N + 1));
    for (m = 0; m < N; m++)
        pts[m] = (fmpz_cmp_ui(smk_points + m, J) > 0) ? J + 1 : fmpz_get_ui(smk_points + m);
    pts[N] = UWORD_MAX;

    a.pts = pts;
    a.A = A;
    a.B = B;
    a.K = K;
    a.N = N;
    a.pn = pn;
    a.prec = 128;
    a.J = J;
    a.t0 = t0;
    /* the odd primes up to sqrt(J) */
    {
        ulong sq = n_sqrt(J) + 1, q, x;
        char * comp = flint_calloc(sq + 1, 1);
        small = flint_malloc(sizeof(ulong) * (sq / 2 + 2));
        nsmall = 0;
        for (q = 3; q <= sq; q += 2)
            if (!comp[q])
            {
                small[nsmall++] = q;
                for (x = q * q; x <= sq; x += 2 * q)
                    comp[x] = 1;
            }
        flint_free(comp);
        a.small = small;
        a.nsmall = nsmall;
        /* the table bound: the composites j whose cofactor j/d is at
           most M are products of table entries */
        if (M != 0)
            a.M = FLINT_MIN(J / 2, M);
        else
            a.M = FLINT_MIN(J / 2, (ulong) (PSMK_TABLE_BYTES / (3 * sizeof(d3b_struct))));
        a.M = FLINT_MAX(a.M, sq);
    }
    GR_MUST_SUCCEED(gr_ctx_init_dfloat(Tt, 3, DFLOAT_BALL));
    a.tre = gr_heap_init_vec(a.M + 1, Tt);
    a.tim = gr_heap_init_vec(a.M + 1, Tt);
    a.tlg = gr_heap_init_vec(a.M + 1, Tt);

    /* the chunks: the dyadic ranges [2^i, 2^(i+1)) while 2^i <= M (the
       cofactors of a range, at most j/2, lie in the previous ones),
       then [2^i, J]; each range in chunks of a fixed length */
    {
        ulong clen = FLINT_MAX(UWORD(1) << 18, J / 4096 + 1), lo, hi, x;
        slong cap = 64 + 2 * (J / clen + 1);
        a.cstart = cs = flint_malloc(sizeof(ulong) * cap);
        a.cend = ce = flint_malloc(sizeof(ulong) * cap);
        rstart = flint_malloc(sizeof(slong) * 66);
        nchunks = 0;
        nranges = 0;
        for (lo = 1; lo <= J; lo = hi + 1)
        {
            hi = (lo <= a.M && lo <= J / 2) ? 2 * lo - 1 : J;
            hi = FLINT_MIN(hi, J);
            rstart[nranges++] = nchunks;
            for (x = lo; x <= hi; x += clen)
            {
                cs[nchunks] = x;
                ce[nchunks] = FLINT_MIN(hi, x + clen - 1);
                nchunks++;
            }
        }
        rstart[nranges] = nchunks;
    }
    a.m0 = flint_malloc(sizeof(slong) * nchunks);
    a.m1 = flint_malloc(sizeof(slong) * nchunks);
    a.res = flint_calloc(nchunks, sizeof(acb_ptr));
    a.status = flint_calloc(nchunks, sizeof(int));

    /* the ranges in order; each in waves of a few chunks per thread,
       whose results are added into S5 (in the order of the chunks)
       before the next wave */
    for (m = 0; m < K * N * 5; m++)
        S5[m] = 0.0;
    acb_init(y);
    {
        slong Wn = 2 * FLINT_MAX(1, flint_get_num_threads()), cw, nw;
        for (r = 0; r < nranges; r++)
        {
            for (cw = rstart[r]; cw < rstart[r + 1]; cw += nw)
            {
                nw = FLINT_MIN(Wn, rstart[r + 1] - cw);
                a.c0 = cw;
                flint_parallel_do(psmk_chunk, &a, nw, 0, FLINT_PARALLEL_STRIDED);
                for (c = cw; c < cw + nw; c++)
                {
                    slong nb = a.m1[c] - a.m0[c] + 1;
                    status |= a.status[c];
                    if (status == GR_SUCCESS)
                    {
                        for (m = a.m0[c]; m <= a.m1[c]; m++)
                        {
                            for (k = 0; k < K; k++)
                            {
                                double * e = S5 + (k * N + m) * 5;
                                psmk_acb_from5(y, e);
                                acb_add(y, y, a.res[c] + (m - a.m0[c]) * K + k, 256);
                                psmk_acb_to5(e, y);
                            }
                        }
                    }
                    _acb_vec_clear(a.res[c], nb * K);
                    a.res[c] = NULL;
                }
            }
        }
    }
    acb_clear(y);

    flint_free(pts);
    flint_free(a.m0);
    flint_free(a.m1);
    flint_free(a.res);
    flint_free(a.status);
    flint_free(small);
    flint_free(cs);
    flint_free(ce);
    flint_free(rstart);
    gr_heap_clear_vec(a.tre, a.M + 1, Tt);
    gr_heap_clear_vec(a.tim, a.M + 1, Tt);
    gr_heap_clear_vec(a.tlg, a.M + 1, Tt);
    gr_ctx_clear(Tt);
    return status == GR_SUCCESS;
}

int
_dfloat_platt_smk(acb_ptr S_table, const fmpz * smk_points, const arb_t t0,
    slong A, slong B, ulong J, slong K, int pn, ulong M, slong prec)
{
    slong N = A * B, i;
    double * S5;
    acb_t cc;
    arb_t w;
    int ok;

    if (pn < 3 || pn > 4 || !dfloat_is_supported() || J < 1)
        return 0;

    S5 = flint_malloc(sizeof(double) * 5 * K * N);
    ok = _dfloat_platt_smk_dd(S5, smk_points, t0, A, B, J, K, pn, M);
    if (ok)
    {
        /* S_table = c S, c = exp(-i t0 log sqrt(pi)) */
        acb_init(cc);
        arb_init(w);
        arb_const_sqrt_pi(w, prec + 64);
        arb_log(w, w, prec + 64);
        arb_mul(w, w, t0, prec + 64);
        arb_neg(w, w);
        arb_sin_cos(acb_imagref(cc), acb_realref(cc), w, prec);
        for (i = 0; i < K * N; i++)
        {
            psmk_acb_from5(S_table + i, S5 + 5 * i);
            acb_mul(S_table + i, S_table + i, cc, prec);
        }
        acb_clear(cc);
        arb_clear(w);
    }
    flint_free(S5);
    return ok;
}
