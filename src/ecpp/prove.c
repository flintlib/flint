/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "thread_support.h"
#include "thread_pool.h"
#include "acb_modular.h"
#include "ecpp.h"

/*
    Atkin-Morain ECPP with the fastECPP organisation of the discriminant
    search (Franke, Kleinjung, Morain, Wirth 2004; Morain 2007), in the
    form of the implementations in PARI/GP (primecert, by Jared Asuncion)
    and in Andreas Enge's CM library (see the documentation for what was
    taken from each; the parameters were tuned against both).

    For each n in the chain:
      1. run through the pool of fundamental discriminants D that are
         smooth over a fixed set of small primes (see disc.c), skipping
         those with a genus character (p^* / n) = -1 (Cornacchia cannot
         succeed for them); sqrt(D) mod n is the product of the square
         roots of the p^*, each computed once per n; Cornacchia writes
         4n = t^2 + |D| v^2 and yields the possible curve orders;
      2. once a batch of orders is collected, remove all their prime
         factors below 2^tdivexp at once with a remainder tree, keep the
         cofactors q with (n^{1/4} + 1)^2 < q <= n/2, sort them by size
         and test them for probable primality in that order: the first
         probable prime is the smallest, i.e. the largest gain in bits
         that the batch offers;
      3. realise the candidate: a root of the Hilbert class polynomial
         H_D modulo n gives the curve, the right twist and a point P
         with m P = O and (m/q) P != O are found by testing; recurse on q.
         A dead end continues with the next probable prime of the batch,
         then with the next batch.
*/

#ifdef ECPP_PROFILE
#include <time.h>
static double _now(void)
{
    struct timespec t;
    timespec_get(&t, TIME_UTC);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}
static double prof_scan, prof_hilbert, prof_roots, prof_points;
static double prof_sqrt, prof_corn, prof_split, prof_bpsw;
static slong prof_hsum, prof_nsteps, prof_twists, prof_gain;
static slong prof_nsqrt, prof_ncorn, prof_ncornok, prof_nsplit, prof_nbpsw, prof_nbatch, prof_ndisc;
static slong prof_nobatch, prof_nobatch_tests, prof_fullbatch;
#define PROF_START(v) double v = _now()
#define PROF_ADD(acc, v) acc += _now() - v
#else
#define PROF_START(v)
#define PROF_ADD(acc, v)
#endif

typedef struct
{
    slong D;
    slong h;
    slong o;                            /* degree of the genus factor of H_D */
    double cost;                        /* estimated realisation cost */
    int use_tower;
    slong Dtower;                       /* D or 4D */
    int veven;                          /* v even in 4n = t^2 + |D| v^2 */
    slong g;
    slong pstar[ECPP_DISC_MAXFAC + 1];  /* the p_i^* and q0 */
    fmpz sqrts[ECPP_DISC_MAXFAC + 1];   /* their square roots mod n */
    fmpz_t m;
    fmpz_t q;
}
ecpp_cand_struct;

/* runtime verbosity (ecpp_set_verbose) */
int ecpp_verbose = 0;

void
ecpp_set_verbose(int verbose)
{
    ecpp_verbose = verbose;
}

typedef struct
{
    ecpp_disc_struct * discs;
    slong ndiscs;
    ulong * primes;         /* odd primes up to pmax */
    slong nprimes;
    fmpz primorial[32];     /* product of the primes below 2^i, computed on demand */
    int have_primorial[32];
    fmpz primorial_chunk[32][ECPP_PRIMORIAL_CHUNKS]; /* as a product of chunks of
                               equal size, for a parallel reduction */
    flint_rand_t state;
    slong max_tries;        /* candidates realised per n before giving up */
    slong tdivexp;          /* trial division exponent, fixed from the size of
                               the input for the whole chain (as in CM): the
                               gain per step does not decrease down the
                               chain, and the primorial is computed once */
    /*
        The realisation of a step (class polynomial, root, curve, point)
        is independent of the search for the next step, which only needs
        the cofactor q: with threads it runs on a pool thread while the
        search for q proceeds. At most one realisation is pending; it is
        joined by the next level before that level starts its own.
    */
    struct
    {
        int active;             /* a realisation is pending */
        int threaded;           /* on a pool thread (else already done) */
        slong cert_index;       /* the step's position in the certificate */
        thread_pool_handle handle;
        ecpp_step_struct step;
        ecpp_cand_struct cand;
        fmpz_t n;
        fmpz_mod_ctx_struct * mctx;
        flint_rand_t state;
        int result;
    }
    pending;
    int * results;          /* result of the realisation, per certificate index */
    slong results_alloc;
    int use_threads;        /* threads available and n large enough */
}
ecpp_ctx_struct;



/*
    Size dependent parameters: bits, largest prime of the set, class number
    bound, trial division exponent. The class number bound and the trial
    division exponent follow PARI/GP (with one more bit of trial division
    from 2000 bits on: the batch factoring gets twice as expensive but the
    batches are twice as large and the gain per step grows by a few bits,
    which was measured to pay). The prime set is much smaller than PARI's:
    since only discriminants with many small prime factors have a small
    odd class number part, discriminants with a large prime factor are
    rarely useful, while each such prime costs a square root per step;
    at 4000 bits, primes up to 400 instead of 1400 cut the square roots
    per step from 54 to 29 and the total time by a quarter.
*/
static const slong ecpp_tune[][4] =
{
    {  300,  100,   4, 16 },
    {  400,  120,   5, 17 },
    {  600,  150,   6, 18 },
    {  700,  250,   8, 20 },
    { 1000,  350,  14, 21 },
    { 1200,  400,  16, 21 },
    { 1500,  450,  18, 22 },
    { 2000,  450,  22, 23 },
    { 2500,  500,  26, 24 },
    { 3000,  500,  32, 24 },
    { 3500,  500,  38, 25 },
    { 4000,  500,  44, 25 },
    { 4500,  550,  50, 25 },
    { 5000,  600,  56, 26 },
    { 6000,  700,  68, 27 },
    { 7000,  800,  80, 27 },
    { 8000,  900,  92, 28 },
    { 9000, 1000, 104, 29 },
    {10000, 1100, 116, 30 },
};

static void
_ecpp_tune(slong bits, slong * pmax, slong * hmax, slong * tdivexp)
{
    slong i, n = sizeof(ecpp_tune) / sizeof(ecpp_tune[0]);
    for (i = 0; i < n; i++)
        if (bits <= ecpp_tune[i][0])
            break;
    if (i == n)
    {
        i = n - 1;
        *pmax = ecpp_tune[i][1] + (bits - ecpp_tune[i][0]) * 2 / 5;
        *hmax = ecpp_tune[i][2] + (bits - ecpp_tune[i][0]) / 80;
        *tdivexp = FLINT_MIN(30, ecpp_tune[i][3] + (bits - ecpp_tune[i][0]) / 4000);
    }
    else
    {
        *pmax = ecpp_tune[i][1];
        *hmax = ecpp_tune[i][2];
        *tdivexp = ecpp_tune[i][3];
    }
}



/*
    The product of the primes below 2^e. With threads it is also kept as
    a product of ECPP_PRIMORIAL_CHUNKS chunks of nearly equal size (the
    primes in ranges of equal length), for a parallel reduction; the
    chunks are products of the primes of their range, so that the whole
    costs about one primorial.
*/
static const fmpz *
_ecpp_primorial(ecpp_ctx_struct * ctx, slong e)
{
    if (!ctx->have_primorial[e])
    {
        fmpz_init(ctx->primorial + e);
        if (!ctx->use_threads)
            fmpz_primorial(ctx->primorial + e, UWORD(1) << e);
        else
        {
            slong k, np, alloc = 1024;
            ulong p, B = UWORD(1) << e, lim;
            fmpz * v = _fmpz_vec_init(alloc);
            n_primes_t it;

            n_primes_init(it);
            p = n_primes_next(it);
            fmpz_one(ctx->primorial + e);
            for (k = 0; k < ECPP_PRIMORIAL_CHUNKS; k++)
            {
                lim = (k == ECPP_PRIMORIAL_CHUNKS - 1) ? B : (B / ECPP_PRIMORIAL_CHUNKS) * (k + 1);
                np = 0;
                while (p < lim)
                {
                    if (np == alloc)
                    {
                        slong q;
                        alloc *= 2;
                        v = flint_realloc(v, alloc * sizeof(fmpz));
                        for (q = np; q < alloc; q++)
                            v[q] = 0;
                    }
                    fmpz_set_ui(v + np, p);
                    np++;
                    p = n_primes_next(it);
                }
                fmpz_init(ctx->primorial_chunk[e] + k);
                if (np > 0)
                    _fmpz_vec_prod(ctx->primorial_chunk[e] + k, v, np);
                else
                    fmpz_one(ctx->primorial_chunk[e] + k);
                fmpz_mul(ctx->primorial + e, ctx->primorial + e, ctx->primorial_chunk[e] + k);
            }
            n_primes_clear(it);
            _fmpz_vec_clear(v, alloc);
        }
        ctx->have_primorial[e] = 1;
    }
    return ctx->primorial + e;
}

/* r_k = chunk_k mod prod, in parallel */
typedef struct
{
    const fmpz * chunks;
    const fmpz * prod;
    const fmpz_preinvn_struct * inv;
    fmpz * r;
}
_chunk_args;

static void
_chunk_worker(slong k, void * varg)
{
    _chunk_args * args = (_chunk_args *) varg;
    fmpz_fdiv_r_preinvn(args->r + k, args->chunks + k, args->prod, args->inv);
}

/*
    Per-n state of the discriminant search: genus characters and square
    roots of the p^*, position in the pool, and the current batch of
    cofactors sorted by size with the position of the next one to test.
*/
typedef struct
{
    const fmpz * n;
    const fmpz_mod_ctx_struct * mctx;
    slong bits;
    signed char * sym;      /* (p^* / n): 0 unknown, 1, -1 */
    signed char sym2[3];    /* (-4 / n), (8 / n), (-8 / n) */
    fmpz * sqrt;            /* sqrt(p^*) mod n */
    signed char * have;     /* 0 unknown, 1 known, 2 failed (n composite) */
    fmpz_t sqrt2[3];        /* sqrt(-1), sqrt(2), sqrt(-2) */
    signed char have2[3];
    slong pos;              /* next discriminant of the pool (persevere mode) */
    slong pprime;           /* next prime of the set to consider */
    signed char * active;   /* primes with square root computed (or known non-residues) */
    signed char * used;     /* discriminants already tried for this n */
    slong delta;            /* minimum gain in bits */
    fmpz_t ts_r, ts_g;      /* n - 1 = 2^ts_e ts_r, ts_g = z^ts_r of order 2^ts_e
                               for a non-residue z: shared Tonelli-Shanks data */
    slong ts_e;
    ulong nmodp[4];         /* n mod 5, 7, 11, 13: cheap Kummer descents */
    fmpz_t qmin, qmax, np1;
    /* batch: orders m, cofactors q sorted by size, pool index of D, status
       (0 untested, 1 probable prime not yet used, 2 composite or used) */
    fmpz * bm, * bq;
    slong * bD;
    signed char * bstat;
    signed char * bveven;
    slong nb, balloc, bpos;
    /* q's already tried for this n */
    fmpz * tried;
    slong ntried;
    int composite;
}
ecpp_scan_struct;

static void
scan_init(ecpp_scan_struct * s, const fmpz_t n, const fmpz_mod_ctx_t mctx,
                                                        const ecpp_ctx_struct * ctx)
{
    ulong n8;
    slong i;

    s->n = n;
    s->mctx = mctx;
    s->bits = fmpz_bits(n);
    s->sym = flint_calloc(ctx->nprimes, sizeof(signed char));
    s->sqrt = _fmpz_vec_init(ctx->nprimes);
    s->have = flint_calloc(ctx->nprimes, sizeof(signed char));
    for (i = 0; i < 3; i++)
    {
        fmpz_init(s->sqrt2[i]);
        s->have2[i] = 0;
    }
    n8 = fmpz_fdiv_ui(n, 8);
    s->sym2[0] = (n8 % 4 == 1) ? 1 : -1;            /* (-4 / n) */
    s->sym2[1] = (n8 == 1 || n8 == 7) ? 1 : -1;     /* (8 / n)  */
    s->sym2[2] = (n8 == 1 || n8 == 3) ? 1 : -1;     /* (-8 / n) */
    s->pos = 0;
    s->pprime = 0;
    s->nmodp[0] = fmpz_fdiv_ui(n, 5);
    s->nmodp[1] = fmpz_fdiv_ui(n, 7);
    s->nmodp[2] = fmpz_fdiv_ui(n, 11);
    s->nmodp[3] = fmpz_fdiv_ui(n, 13);
    /* shared Tonelli-Shanks data: one exponentiation per n instead of one
       per square root (and no search for a non-residue per root) */
    fmpz_init(s->ts_r);
    fmpz_init(s->ts_g);
    fmpz_sub_ui(s->ts_r, n, 1);
    s->ts_e = fmpz_val2(s->ts_r);
    fmpz_fdiv_q_2exp(s->ts_r, s->ts_r, s->ts_e);
    {
        ulong z = 2;
        fmpz_t t;
        fmpz_init(t);
        while (1)
        {
            fmpz_set_ui(t, z);
            if (fmpz_jacobi(t, n) == -1)
                break;
            z++;
        }
        fmpz_powm(s->ts_g, t, s->ts_r, n);
        fmpz_clear(t);
    }
    s->active = flint_calloc(ctx->nprimes, sizeof(signed char));
    s->used = flint_calloc(ctx->ndiscs, sizeof(signed char));
    fmpz_init(s->qmin);
    fmpz_init(s->qmax);
    fmpz_init(s->np1);
    /* q must exceed (n^{1/4} + 1)^2; qmin = (floor(n^{1/4}) + 2)^2 suffices */
    fmpz_root(s->qmin, n, 4);
    fmpz_add_ui(s->qmin, s->qmin, 2);
    fmpz_mul(s->qmin, s->qmin, s->qmin);
    /* minimum gain: half the bits removed on average by the trial
       division (Franke, Kleinjung, Morain, Wirth), for large n */
    s->delta = (s->bits >= 1000) ? ctx->tdivexp / 2 + 1 : 1;
    fmpz_fdiv_q_2exp(s->qmax, n, s->delta);
    fmpz_add_ui(s->np1, n, 1);
    s->balloc = 64;
    s->bm = _fmpz_vec_init(s->balloc);
    s->bq = _fmpz_vec_init(s->balloc);
    s->bD = flint_malloc(s->balloc * sizeof(slong));
    s->bstat = flint_malloc(s->balloc * sizeof(signed char));
    s->bveven = flint_malloc(s->balloc * sizeof(signed char));
    s->nb = s->bpos = 0;
    s->tried = _fmpz_vec_init(ctx->max_tries + 2);
    s->ntried = 0;
    s->composite = 0;
}

static void
scan_clear(ecpp_scan_struct * s, const ecpp_ctx_struct * ctx)
{
    slong i;
    flint_free(s->sym);
    _fmpz_vec_clear(s->sqrt, ctx->nprimes);
    flint_free(s->have);
    for (i = 0; i < 3; i++)
        fmpz_clear(s->sqrt2[i]);
    flint_free(s->active);
    flint_free(s->used);
    fmpz_clear(s->ts_r);
    fmpz_clear(s->ts_g);
    fmpz_clear(s->qmin);
    fmpz_clear(s->qmax);
    fmpz_clear(s->np1);
    _fmpz_vec_clear(s->bm, s->balloc);
    _fmpz_vec_clear(s->bq, s->balloc);
    flint_free(s->bD);
    flint_free(s->bstat);
    flint_free(s->bveven);
    _fmpz_vec_clear(s->tried, ctx->max_tries + 2);
}

static void
scan_batch_fit(ecpp_scan_struct * s, slong len)
{
    if (len > s->balloc)
    {
        slong i, old = s->balloc;
        s->balloc = FLINT_MAX(len, 2 * old);
        s->bm = flint_realloc(s->bm, s->balloc * sizeof(fmpz));
        s->bq = flint_realloc(s->bq, s->balloc * sizeof(fmpz));
        s->bD = flint_realloc(s->bD, s->balloc * sizeof(slong));
        s->bstat = flint_realloc(s->bstat, s->balloc * sizeof(signed char));
        s->bveven = flint_realloc(s->bveven, s->balloc * sizeof(signed char));
        for (i = old; i < s->balloc; i++)
            s->bm[i] = s->bq[i] = 0;    /* fresh (uninitialised) memory: cannot use fmpz_zero */
    }
}

/* parallel computation of the missing square roots for a set of discriminants */
typedef struct
{
    ecpp_scan_struct * s;
    const ulong * primes;
    slong * which;      /* prime indices; -1, -2, -3 stand for -1, 2, -2 */
}
_sqrt_args;

/* square root of a (a square) modulo n with the shared Tonelli-Shanks data */
static int
_sqrtmod_shared(fmpz_t x, const fmpz_t a, const ecpp_scan_struct * s)
{
    const fmpz * n = s->n;
    fmpz_t b, g, t, e;
    slong m, r, i;
    int ok = 1;

    fmpz_init(b); fmpz_init(g); fmpz_init(t); fmpz_init(e);
    /* x = a^{(r+1)/2}, b = a^r = x^2 / a */
    fmpz_add_ui(e, s->ts_r, 1);
    fmpz_fdiv_q_2exp(e, e, 1);
    fmpz_powm(x, a, e, n);
    fmpz_mod_mul(b, x, x, s->mctx);
    fmpz_mod_inv(t, a, s->mctx);
    fmpz_mod_mul(b, b, t, s->mctx);
    fmpz_set(g, s->ts_g);
    r = s->ts_e;
    while (!fmpz_is_one(b))
    {
        /* smallest m with b^(2^m) = 1 */
        fmpz_set(t, b);
        for (m = 0; m < r && !fmpz_is_one(t); m++)
            fmpz_mod_mul(t, t, t, s->mctx);
        if (m == r)
        {
            ok = 0;
            break;
        }
        /* t = g^(2^(r-m-1)) */
        fmpz_set(t, g);
        for (i = 0; i < r - m - 1; i++)
            fmpz_mod_mul(t, t, t, s->mctx);
        fmpz_mod_mul(g, t, t, s->mctx);
        fmpz_mod_mul(x, x, t, s->mctx);
        fmpz_mod_mul(b, b, g, s->mctx);
        r = m;
    }
    fmpz_clear(b); fmpz_clear(g); fmpz_clear(t); fmpz_clear(e);
    return ok;
}

static void
_sqrt_worker(slong i, void * varg)
{
    _sqrt_args * args = (_sqrt_args *) varg;
    ecpp_scan_struct * s = args->s;
    slong k = args->which[i];
    fmpz_t t;
    int ok;

    fmpz_init(t);
    if (k < 0)
    {
        if (k == -1) fmpz_sub_ui(t, s->n, 1);
        else if (k == -2) fmpz_set_ui(t, 2);
        else fmpz_sub_ui(t, s->n, 2);
        ok = _sqrtmod_shared(s->sqrt2[-k - 1], t, s);
        s->have2[-k - 1] = ok ? 1 : 2;
    }
    else
    {
        ulong p = args->primes[k];
        if (p % 4 == 1)
            fmpz_set_ui(t, p);
        else
            fmpz_sub_ui(t, s->n, p);
        ok = _sqrtmod_shared(s->sqrt + k, t, s);
        s->have[k] = ok ? 1 : 2;
    }
    fmpz_clear(t);
}

static void
_scan_sqrt_prefetch(ecpp_scan_struct * s, const ecpp_ctx_struct * ctx,
                                                const slong * idx, slong num)
{
    slong i, j, count = 0;
    slong * which = flint_malloc((ctx->nprimes + 3) * sizeof(slong));
    signed char * need = flint_calloc(ctx->nprimes + 3, sizeof(signed char));
    _sqrt_args args;

    for (i = 0; i < num; i++)
    {
        const ecpp_disc_struct * d = ctx->discs + idx[i];
        if (d->q0 == -4) need[ctx->nprimes + 0] = 1;
        if (d->q0 == 8) need[ctx->nprimes + 1] = 1;
        if (d->q0 == -8) need[ctx->nprimes + 2] = 1;
        for (j = 0; j < d->nfac; j++)
            need[d->fac[j]] = 1;
    }
    for (j = 0; j < 3; j++)
        if (need[ctx->nprimes + j] && s->have2[j] == 0)
            which[count++] = -1 - j;
    for (j = 0; j < ctx->nprimes; j++)
        if (need[j] && s->have[j] == 0)
            which[count++] = j;

    args.s = s;
    args.primes = ctx->primes;
    args.which = which;
    if (count == 1 || !ctx->use_threads)
    {
        slong q;
        for (q = 0; q < count; q++)
            _sqrt_worker(q, &args);
    }
    else if (count > 1)
        flint_parallel_do(_sqrt_worker, &args, count, 0, FLINT_PARALLEL_STRIDED);
#ifdef ECPP_PROFILE
    prof_nsqrt += count;
#endif

    flint_free(which);
    flint_free(need);
}

/* sqrt(D) mod n from the cached roots; 0 if some root is missing */
static int
_scan_sqrt_disc(fmpz_t r, const ecpp_scan_struct * s, const ecpp_disc_struct * d)
{
    slong i;

    if (d->q0 == 1)
        fmpz_one(r);
    else
    {
        slong j = (d->q0 == -4) ? 0 : (d->q0 == 8) ? 1 : 2;
        if (s->have2[j] != 1)
            return 0;
        fmpz_mod_add(r, s->sqrt2[j], s->sqrt2[j], s->mctx);
    }
    for (i = 0; i < d->nfac; i++)
    {
        slong k = d->fac[i];
        if (s->have[k] != 1)
            return 0;
        fmpz_mod_mul(r, r, s->sqrt + k, s->mctx);
    }
    return 1;
}

/* Cornacchia and the possible orders for a set of discriminants, in parallel */
typedef struct
{
    ecpp_scan_struct * s;
    const ecpp_ctx_struct * ctx;
    const slong * idx;
    fmpz * m;           /* 6 per discriminant */
    slong * norders;
    signed char * veven;
}
_corn_args;

static void
_corn_worker(slong w, void * varg)
{
    _corn_args * args = (_corn_args *) varg;
    ecpp_scan_struct * s = args->s;
    const ecpp_disc_struct * d = args->ctx->discs + args->idx[w];
    fmpz * ms = args->m + 6 * w;
    fmpz_t sqrtD, t, v, u;
    slong norders = 0;

    fmpz_init(sqrtD); fmpz_init(t); fmpz_init(v); fmpz_init(u);

    if (_scan_sqrt_disc(sqrtD, s, d) && ecpp_cornacchia(t, v, s->n, d->D, sqrtD))
    {
        fmpz_add((ms + 0), s->np1, t);
        fmpz_sub((ms + 1), s->np1, t);
        norders = 2;
        if (d->D == -4)
        {
            fmpz_mul_2exp(u, v, 1);
            fmpz_add((ms + 2), s->np1, u);
            fmpz_sub((ms + 3), s->np1, u);
            norders = 4;
        }
        else if (d->D == -3)
        {
            fmpz_mul_ui(u, v, 3);
            fmpz_add((ms + 2), t, u);
            fmpz_fdiv_q_2exp((ms + 2), (ms + 2), 1);
            fmpz_sub((ms + 3), t, u);
            fmpz_fdiv_q_2exp((ms + 3), (ms + 3), 1);
            fmpz_add((ms + 4), s->np1, (ms + 2));
            fmpz_sub((ms + 5), s->np1, (ms + 2));
            fmpz_add((ms + 2), s->np1, (ms + 3));
            fmpz_sub((ms + 3), s->np1, (ms + 3));
            norders = 6;
        }
    }
    args->norders[w] = norders;
    args->veven[w] = (norders > 0) ? fmpz_is_even(v) : 0;

    fmpz_clear(sqrtD); fmpz_clear(t); fmpz_clear(v); fmpz_clear(u);
}

typedef struct
{
    fmpz q;
    fmpz m;
    slong D;
    int veven;
}
_bq_entry;

static int
_bq_cmp(const void * x, const void * y)
{
    return fmpz_cmp(&((const _bq_entry *) x)->q, &((const _bq_entry *) y)->q);
}

/*
    Collects the next batch of orders. As in Enge's CM library (following
    Franke, Kleinjung, Morain and Wirth), the set of primes with a square
    root modulo n is extended, in increasing order and skipping the
    non-residues, until the discriminants of the pool that become
    available (all their odd primes in the set, the character at 2 equal
    to 1, not yet used for this n) are expected to yield about three
    prime cofactors; then square roots (in parallel), Cornacchia for all
    those discriminants (in parallel), removal of the small prime factors
    of the orders with one remainder tree, and the admissible cofactors
    sorted by size. Returns the number of cofactors, 0 when the primes are
    exhausted.
*/
static slong
_scan_fill_batch(ecpp_scan_struct * s, ecpp_ctx_struct * ctx)
{
    slong pmax_tune, hmax, tdivexp, i, j, w;
    slong * idx, * norders;
    fmpz * ms;
    signed char * veven, * mveven;
    _corn_args cargs;
    slong nm = 0, malloc_m = 64, nidx = 0, idx_alloc = 256;
    fmpz * mlist;
    slong * mD;
    double prob_prime, expected;
/* expected prime cofactors per round; the estimate below is optimistic
   (the admissibility of a cofactor is not accounted for): with 3.0 about
   11% of the rounds yield no prime and cost a full batch of tests */
    const double min_prime = ECPP_MIN_PRIME;

    _ecpp_tune(s->bits, &pmax_tune, &hmax, &tdivexp);
    tdivexp = ctx->tdivexp;
    if (s->bits < 1000)
    {
        /* down the chain of a small input, the primorial's size matters:
           the size-tuned bound rather than the input's */
        slong pm2, hm2;
        _ecpp_tune(s->bits, &pm2, &hm2, &tdivexp);
        tdivexp = FLINT_MIN(tdivexp, ctx->tdivexp);
    }
    prob_prime = 1.7811 * tdivexp / (double) s->bits;

    idx = flint_malloc(idx_alloc * sizeof(slong));
    mlist = _fmpz_vec_init(malloc_m);
    mD = flint_malloc(malloc_m * sizeof(slong));
    mveven = flint_malloc(malloc_m * sizeof(signed char));

    s->nb = s->bpos = 0;

    while (nm == 0 && s->pprime < ctx->nprimes)
    {
        slong * newp = flint_malloc(ctx->nprimes * sizeof(slong));
        slong nnew = 0;

        /* extend the prime set until enough is expected */
        expected = 0.0;
        while (expected < min_prime && s->pprime < ctx->nprimes)
        {
            slong k = s->pprime++;
            ulong p = ctx->primes[k];

            if (s->sym[k] == 0)
                s->sym[k] = n_jacobi(fmpz_fdiv_ui(s->n, p), p) == 1 ? 1 : -1;
            if (s->sym[k] != 1)
                continue;
            s->active[k] = 1;
            newp[nnew++] = k;

            /* discriminants that become available with p */
            for (i = 0; i < ctx->ndiscs; i++)
            {
                const ecpp_disc_struct * d = ctx->discs + i;
                int ok = !s->used[i], has = 0;
                if (d->q0 == -4 && s->sym2[0] != 1) ok = 0;
                if (d->q0 == 8 && s->sym2[1] != 1) ok = 0;
                if (d->q0 == -8 && s->sym2[2] != 1) ok = 0;
                for (j = 0; j < d->nfac && ok; j++)
                {
                    if (d->fac[j] == k)
                        has = 1;
                    else if (!s->active[d->fac[j]])
                        ok = 0;
                }
                if (!ok || !has)
                    continue;
                s->used[i] = 1;
                if (nidx == idx_alloc)
                {
                    idx_alloc *= 2;
                    idx = flint_realloc(idx, idx_alloc * sizeof(slong));
                }
                idx[nidx++] = i;
                /* orders per success times the success probability 2^{g-1}/h */
                expected += prob_prime * ((d->D == -3) ? 6.0 : (d->D == -4) ? 4.0 : 2.0) / d->o;
            }
        }
        /* discriminants without odd prime factor (-4, -8) on the first round */
        if (s->pprime <= nnew + 1)
        {
            for (i = 0; i < ctx->ndiscs; i++)
            {
                const ecpp_disc_struct * d = ctx->discs + i;
                if (d->nfac == 0 && !s->used[i])
                {
                    if ((d->q0 == -4 && s->sym2[0] == 1) || (d->q0 == 8 && s->sym2[1] == 1)
                                                        || (d->q0 == -8 && s->sym2[2] == 1))
                    {
                        s->used[i] = 1;
                        idx[nidx++] = i;
                    }
                }
            }
        }
        flint_free(newp);
        if (nidx == 0)
            continue;
#ifdef ECPP_PROFILE
        prof_ndisc += nidx;
#endif

        {
            PROF_START(ta);
            _scan_sqrt_prefetch(s, ctx, idx, nidx);
            PROF_ADD(prof_sqrt, ta);
        }

        norders = flint_malloc(nidx * sizeof(slong));
        ms = _fmpz_vec_init(6 * nidx);
        veven = flint_malloc(nidx * sizeof(signed char));
        {
            PROF_START(tb);
            cargs.s = s; cargs.ctx = ctx; cargs.idx = idx; cargs.m = ms; cargs.norders = norders;
            cargs.veven = veven;
            if (nidx == 1 || !ctx->use_threads)
            {
                for (w = 0; w < nidx; w++)
                    _corn_worker(w, &cargs);
            }
            else
                flint_parallel_do(_corn_worker, &cargs, nidx, 0, FLINT_PARALLEL_STRIDED);
            PROF_ADD(prof_corn, tb);
        }

        for (w = 0; w < nidx; w++)
        {
            for (j = 0; j < norders[w]; j++)
            {
                if (fmpz_sgn(ms + 6 * w + j) <= 0)
                    continue;
                if (nm == malloc_m)
                {
                    malloc_m *= 2;
                    mlist = flint_realloc(mlist, malloc_m * sizeof(fmpz));
                    mD = flint_realloc(mD, malloc_m * sizeof(slong));
                    mveven = flint_realloc(mveven, malloc_m * sizeof(signed char));
                    for (i = nm; i < malloc_m; i++)
                        mlist[i] = 0;
                }
                fmpz_set(mlist + nm, ms + 6 * w + j);
                mD[nm] = idx[w];
                mveven[nm] = veven[w];
                nm++;
            }
        }
        flint_free(norders);
        flint_free(veven);
        _fmpz_vec_clear(ms, 6 * nidx);
        nidx = 0;
    }

    if (nm > 0)
    {
        /* batch removal of the primes below 2^tdivexp: r_i = primorial mod m_i */
        fmpz_multi_mod_t P;
        fmpz * r;
        fmpz_t g, q;
        _bq_entry * ent;
        slong ne = 0;
        PROF_START(tc);

        r = _fmpz_vec_init(nm);
        /* the primorial is huge: reduce it once modulo the product of the
           orders (fmpz_multi_mod_precomp would do this twice), then the
           remainder tree */
        fmpz_multi_mod_init(P);
        fmpz_multi_mod_precompute(P, mlist, nm);
        {
            /* the big division with a precomputed inverse (about three
               times faster than GMP's division at these sizes) */
            fmpz_t prod, r0;
            fmpz_preinvn_t inv;
            fmpz_init(prod);
            fmpz_init(r0);
            _fmpz_vec_prod(prod, mlist, nm);
            fmpz_preinvn_init(inv, prod);
            {
                /* the primorial as a product of chunks: reduced in
                   parallel, then multiplied modulo prod */
                _chunk_args cargs2;
                fmpz * rk = _fmpz_vec_init(ECPP_PRIMORIAL_CHUNKS);
                slong k;
                _ecpp_primorial(ctx, tdivexp);
                cargs2.chunks = ctx->primorial_chunk[tdivexp];
                cargs2.prod = prod;
                cargs2.inv = inv;
                cargs2.r = rk;
                if (ctx->use_threads)
                {
                    flint_parallel_do(_chunk_worker, &cargs2, ECPP_PRIMORIAL_CHUNKS, 0, FLINT_PARALLEL_STRIDED);
                    fmpz_set(r0, rk + 0);
                    for (k = 1; k < ECPP_PRIMORIAL_CHUNKS; k++)
                    {
                        fmpz_mul(r0, r0, rk + k);
                        fmpz_fdiv_r_preinvn(r0, r0, prod, inv);
                    }
                }
                else
                    fmpz_fdiv_r_preinvn(r0, ctx->primorial + tdivexp, prod, inv);
                _fmpz_vec_clear(rk, ECPP_PRIMORIAL_CHUNKS);
            }
            fmpz_preinvn_clear(inv);
            fmpz_multi_mod_precomp(r, P, r0, 0);
            fmpz_clear(prod);
            fmpz_clear(r0);
        }
        fmpz_multi_mod_clear(P);

        fmpz_init(g);
        fmpz_init(q);
        ent = flint_malloc(nm * sizeof(_bq_entry));

        for (i = 0; i < nm; i++)
        {
            fmpz_gcd(g, r + i, mlist + i);
            if (fmpz_is_one(g))
                continue;       /* no small factor at all: unusable */
            fmpz_set(q, mlist + i);
            while (!fmpz_is_one(g))
            {
                fmpz_divexact(q, q, g);
                fmpz_gcd(g, q, g);
            }
            if (fmpz_cmp(q, s->qmin) < 0 || fmpz_cmp(q, s->qmax) > 0)
                continue;
            for (j = 0; j < s->ntried; j++)
                if (fmpz_equal(q, s->tried + j))
                    break;
            if (j < s->ntried)
                continue;
            fmpz_init_set(&ent[ne].q, q);
            fmpz_init_set(&ent[ne].m, mlist + i);
            ent[ne].D = mD[i];
            ent[ne].veven = mveven[i];
            ne++;
        }
        qsort(ent, ne, sizeof(_bq_entry), _bq_cmp);

        scan_batch_fit(s, ne);
        for (i = 0; i < ne; i++)
        {
            fmpz_swap(s->bq + i, &ent[i].q);
            fmpz_swap(s->bm + i, &ent[i].m);
            s->bD[i] = ent[i].D;
            s->bveven[i] = ent[i].veven;
            s->bstat[i] = 0;
            fmpz_clear(&ent[i].q);
            fmpz_clear(&ent[i].m);
        }
        s->nb = ne;
        flint_free(ent);
        _fmpz_vec_clear(r, nm);
        fmpz_clear(g);
        fmpz_clear(q);
        PROF_ADD(prof_split, tc);
#ifdef ECPP_PROFILE
        prof_nsplit += nm; prof_nbatch++;
#endif
    }

    flint_free(idx);
    _fmpz_vec_clear(mlist, malloc_m);
    flint_free(mD);
    flint_free(mveven);

    return s->nb;
}

/* probable primality of a chunk of cofactors, in parallel */
typedef struct
{
    const fmpz * q;
    int * res;
}
_bpsw_args;

static void
_bpsw_worker(slong i, void * varg)
{
    _bpsw_args * args = (_bpsw_args *) varg;
    /* the cofactors have no small prime factors: straight to BPSW */
    args->res[i] = fmpz_is_probabprime_BPSW(args->q + i);
}

/*
    The next candidate (D, m, q) for n: among the probable prime cofactors
    of the current batch (all of which are tested, in parallel, when the
    batch is new), the one with the best estimated cost per bit gained that
    has not been returned yet; when the batch is used up, the next batch.
    Returns 1 and fills c, or 0 if the pool is exhausted.
*/
static int
_scan_next_candidate(ecpp_scan_struct * s, ecpp_ctx_struct * ctx, ecpp_cand_struct * c)
{
    slong nthreads = ctx->use_threads ? flint_get_num_threads() : 1;
    int * res = flint_malloc(nthreads * sizeof(int));
    _bpsw_args args;
    int found = 0;

    /* tests beyond the first probable prime when its degree is not small */
    const slong extra = FLINT_MAX(16, 4 * nthreads);

    while (!found)
    {
        slong i, best = -1;
        double best_score = 0;

        {
            int unused = 0;
            for (i = 0; i < s->bpos; i++)
                if (s->bstat[i] == 1)
                    unused = 1;
            if (s->bpos >= s->nb && !unused)
            {
                /* current batch used up: next one */
                if (s->pprime >= ctx->nprimes)
                {
                    /*
                        The prime set is exhausted without a candidate
                        (an n for which few discriminants are usable, e.g.
                        n = 7 mod 8 excludes all D = -4m, -8m, and no
                        cofactor happened to be prime). Before giving up,
                        run the pool again with every gain accepted: the
                        minimum gain is only a speed heuristic. The square
                        roots are kept.
                    */
                    if (s->delta > 1)
                    {
                        s->delta = 1;
                        fmpz_fdiv_q_2exp(s->qmax, s->n, 1);
                        s->pprime = 0;
                        s->pos = 0;
                        memset(s->active, 0, ctx->nprimes);
                        memset(s->used, 0, ctx->ndiscs);
                        s->nb = s->bpos = 0;
                        continue;
                    }
                    break;
                }
                if (_scan_fill_batch(s, ctx) == 0)
                    continue;
            }
        }

        /*
            Test the untested cofactors, smallest first. The first probable
            prime is the one with the largest gain; testing further only
            pays if a candidate with a much cheaper realisation turns up,
            so stop as soon as one with a cheap realisation has been
            found, and otherwise after a few more tests.
        */
        while (s->bpos < s->nb)
        {
            if (best >= 0 && (ecpp_disc_cost_n(ctx->discs + s->bD[best], s->bveven[best], s->nmodp, s->bits) < 6.0
                                || s->bpos >= best + extra))
                break;
            slong chunk = FLINT_MIN(nthreads, s->nb - s->bpos), b;
            PROF_START(td);

            args.q = s->bq + s->bpos;
            args.res = res;
            if (chunk == 1 || !ctx->use_threads)
            {
                for (b = 0; b < chunk; b++)
                    _bpsw_worker(b, &args);
            }
            else
                flint_parallel_do(_bpsw_worker, &args, chunk, 0, FLINT_PARALLEL_STRIDED);
            PROF_ADD(prof_bpsw, td);
#ifdef ECPP_PROFILE
            prof_nbpsw += chunk;
#endif
            for (b = 0; b < chunk; b++)
            {
                s->bstat[s->bpos + b] = res[b] ? 1 : 2;
                if (res[b] && best == -1)
                    best = s->bpos + b;
            }
            s->bpos += chunk;
        }

#ifdef ECPP_PROFILE
        {
            slong np = 0;
            for (i = 0; i < s->bpos; i++)
                np += (s->bstat[i] == 1);
            if (np == 0 && s->bpos >= s->nb)
                prof_nobatch++, prof_nobatch_tests += s->nb;
            if (s->bpos == s->nb) prof_fullbatch++;
        }
#endif
        best = -1;
        for (i = 0; i < s->bpos; i++)
        {
            if (s->bstat[i] == 1)
            {
                slong gain = FLINT_MAX(1, s->bits - (slong) fmpz_bits(s->bq + i));
                /* realisation cost plus the per-step cost of the scan
                   (square roots, batch factoring, tests), in units of a
                   scalar multiplication, per bit gained */
                double cst = ecpp_disc_cost_n(ctx->discs + s->bD[i], s->bveven[i], s->nmodp, s->bits);
                double score = (cst + 4.0 + s->bits / 500.0) / gain;
                if (cst >= 1e5)
                {
                    s->bstat[i] = 2;    /* never realised (no way to a curve) */
                    continue;
                }
                if (best == -1 || score < best_score)
                {
                    best = i;
                    best_score = score;
                }
            }
        }

        if (best >= 0)
        {
            const ecpp_disc_struct * d = ctx->discs + s->bD[best];
            c->D = d->D;
            c->h = d->h;
            c->o = d->o;
            c->cost = d->cost;
            c->veven = s->bveven[best];
            c->use_tower = ecpp_disc_use_tower(d, c->veven, &c->Dtower);
            c->cost = ecpp_disc_cost_n(d, c->veven, s->nmodp, s->bits);
            c->g = 0;
            for (i = 0; i < d->nfac; i++)
            {
                ulong p = ctx->primes[d->fac[i]];
                c->pstar[c->g] = (p % 4 == 1) ? (slong) p : -(slong) p;
                fmpz_set(c->sqrts + c->g, s->sqrt + d->fac[i]);
                c->g++;
            }
            if (d->q0 != 1)
            {
                slong j = (d->q0 == -4) ? 0 : (d->q0 == 8) ? 1 : 2;
                c->pstar[c->g] = d->q0;
                fmpz_mod_add(c->sqrts + c->g, s->sqrt2[j], s->sqrt2[j], s->mctx);
                c->g++;
            }
            fmpz_set(c->m, s->bm + best);
            fmpz_set(c->q, s->bq + best);
            s->bstat[best] = 2;
            found = 1;
        }
    }

    flint_free(res);
    return found;
}

/* the curve y^2 = x^3 + a x + b with given j-invariant (j != 0, 1728) */
static void
_ecpp_curve_from_j(fmpz_t a, fmpz_t b, const fmpz_t j, const fmpz_mod_ctx_t ctx)
{
    fmpz_t k, t;

    fmpz_init(k);
    fmpz_init(t);

    /* k = j / (1728 - j); a = 3k, b = 2k */
    fmpz_set_ui(t, 1728);
    fmpz_mod_sub(t, t, j, ctx);
    fmpz_mod_inv(t, t, ctx);
    fmpz_mod_mul(k, j, t, ctx);
    fmpz_mod_add(a, k, k, ctx);
    fmpz_mod_add(b, a, k, ctx);     /* b = 3k for the moment */
    fmpz_swap(a, b);                /* a = 3k, b = 2k */

    fmpz_clear(k);
    fmpz_clear(t);
}

/* a random point on y^2 = x^3 + a x + b; returns 0 if none was found */
static int
_ecpp_random_point(fmpz_t x, fmpz_t y, const fmpz_t a, const fmpz_t b,
                                    flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t rhs, t;
    slong tries;
    int found = 0;

    fmpz_init(rhs);
    fmpz_init(t);

    for (tries = 0; tries < 200 && !found; tries++)
    {
        fmpz_randm(x, state, n);
        fmpz_mod_mul(rhs, x, x, ctx);
        fmpz_mod_add(rhs, rhs, a, ctx);
        fmpz_mod_mul(rhs, rhs, x, ctx);
        fmpz_mod_add(rhs, rhs, b, ctx);

        if (fmpz_jacobi(rhs, n) != 1)
            continue;
        if (!fmpz_sqrtmod(y, rhs, n))
            continue;
        fmpz_mod_mul(t, y, y, ctx);
        found = fmpz_equal(t, rhs);
    }

    fmpz_clear(rhs);
    fmpz_clear(t);
    return found;
}

/*
    Given a curve (a, b) and m = k q, look for a point P with m P = O and
    k P != O. Returns 1 and sets (x, y) on success, 0 if m is not the
    order of the curve (or the point was unlucky), -1 if n was found to be
    composite. m P is computed as q (k P), with k P normalised to affine
    coordinates first, so that the cost is that of a multiplication by q
    plus one by k.
*/
static int
_ecpp_find_point(fmpz_t x, fmpz_t y, const fmpz_t a, const fmpz_t b,
        const fmpz_t q, const fmpz_t k, flint_rand_t state,
        const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    ecpp_point_t P, R;
    fmpz_t acc, g, zi;
    slong tries;
    int result = 0;

    ecpp_point_init(P);
    ecpp_point_init(R);
    fmpz_init(acc);
    fmpz_init(g);
    fmpz_init(zi);

    for (tries = 0; tries < 4; tries++)
    {
        if (!_ecpp_random_point(x, y, a, b, state, ctx))
        {
            result = 0;
            break;
        }

        ecpp_point_set_affine(P, x, y);
        fmpz_one(acc);

        /* R = k P; then q R should be O and R != O */
        ecpp_point_mul(R, P, k, a, acc, ctx);
        if (ecpp_point_is_zero(R))
            continue;   /* unlucky point (or wrong order); try another */

        /* affine R = (X / Z^2, Y / Z^3) */
        fmpz_gcd(g, R->Z, n);
        if (!fmpz_is_one(g))
        {
            result = -1;
            break;
        }
        fmpz_mod_inv(zi, R->Z, ctx);
        fmpz_mod_mul(g, zi, zi, ctx);
        fmpz_mod_mul(R->X, R->X, g, ctx);
        fmpz_mod_mul(g, g, zi, ctx);
        fmpz_mod_mul(R->Y, R->Y, g, ctx);
        fmpz_one(R->Z);

        ecpp_point_mul(P, R, q, a, acc, ctx);

        fmpz_gcd(g, acc, n);
        if (!fmpz_is_one(g))
        {
            result = -1;
            break;
        }

        if (!ecpp_point_is_zero(P))
        {
            result = 0;     /* m is not the order of this twist */
            break;
        }

        result = 1;
        break;
    }

    ecpp_point_clear(P);
    ecpp_point_clear(R);
    fmpz_clear(acc);
    fmpz_clear(g);
    fmpz_clear(zi);

    return result;
}

/*
    One root modulo n of the Hilbert class polynomial H_D, which splits
    into linear factors modulo n when Cornacchia succeeded for D (n splits
    completely in the Hilbert class field). We therefore skip the
    computation of gcd(x^n - x, H_D) and split directly with random
    (x + r)^{(n-1)/2} - 1, keeping a smaller factor each time.
*/
static int
_ecpp_class_poly_root(fmpz_t j, const ecpp_cand_struct * c, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    slong D = c->D;
    fmpz_poly_t H;
    fmpz_mod_poly_t f;
    int found = 0;

    fmpz_poly_init(H);
    fmpz_mod_poly_init(f, ctx);

    {
        PROF_START(t0);
        /* the class field tower first (Weber invariant for even D), then
           the factor of H_D over the genus field (degree h / 2^{g-1}) when
           there is one, else H_D itself */
        /* odd D with v even: the order of conductor 2 admits the Weber
           invariant, and its curves have Frobenius (t + v sqrt D) / 2 in
           Z[sqrt D], so a curve of the right order is among them */
        if (c->use_tower && ecpp_class_poly_tower(j, c->Dtower, state, ctx))
        {
            PROF_ADD(prof_hilbert, t0);
            found = 1;
            goto cleanup;
        }
        if (c->g < 2 || !ecpp_class_poly_genus(f, D, c->pstar, c->g, c->sqrts, ctx))
        {
            if (c->g >= 2 && c->h > 4 * c->o + 8)
            {
                /* the genus factor failed and H_D is much larger than what
                   the candidate was priced at: give up on this candidate */
                PROF_ADD(prof_hilbert, t0);
                goto cleanup;
            }
            acb_modular_hilbert_class_poly(H, D);
            fmpz_mod_poly_set_fmpz_poly(f, H, ctx);
            fmpz_mod_poly_make_monic(f, f, ctx);
        }
        PROF_ADD(prof_hilbert, t0);
    }

    {
        PROF_START(t1);
        /* a root of f: by radicals up to degree 4, else by random splitting */
        found = ecpp_poly_root(j, f, state, ctx);
        PROF_ADD(prof_roots, t1);
    }

cleanup:
    fmpz_mod_poly_clear(f, ctx);
    fmpz_poly_clear(H);

    return found;
}

/* the two twists of a curve with j != 0, 1728 tried in parallel */
typedef struct
{
    const fmpz * a;
    const fmpz * b;
    const fmpz * g;
    const fmpz * q;
    const fmpz * k;
    const fmpz_mod_ctx_struct * ctx;
    fmpz * x;   /* 2 each */
    fmpz * y;
    fmpz * ta;
    fmpz * tb;
    int * res;
    ulong seed;
}
_twist_args;

static void
_twist_worker(slong i, void * varg)
{
    _twist_args * args = (_twist_args *) varg;
    flint_rand_t state;
    fmpz_t t;

    flint_rand_init(state);
    flint_rand_set_seed(state, args->seed + 1000003 * i, 87654321 + i);
    fmpz_init(t);
    if (i == 0)
    {
        fmpz_set(args->ta + 0, args->a);
        fmpz_set(args->tb + 0, args->b);
    }
    else
    {
        /* (a, b) -> (a g^2, b g^3) */
        fmpz_mod_mul(t, args->g, args->g, args->ctx);
        fmpz_mod_mul(args->ta + 1, args->a, t, args->ctx);
        fmpz_mod_mul(t, t, args->g, args->ctx);
        fmpz_mod_mul(args->tb + 1, args->b, t, args->ctx);
    }
    args->res[i] = _ecpp_find_point(args->x + i, args->y + i, args->ta + i,
                                    args->tb + i, args->q, args->k, state, args->ctx);
    fmpz_clear(t);
    flint_rand_clear(state);
}

/*
    Realise a candidate (D, m, q) for n: find curve and point. Returns 1
    and fills the step, 0 if it did not work out, -1 if n is composite.
*/
static int
_ecpp_realise(ecpp_step_struct * s, const fmpz_t n, const ecpp_cand_struct * c,
                                flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    fmpz_t j, a, b, g, k, t, x, y;
    slong twist, ntwists;
    int result = 0;

    fmpz_init(j); fmpz_init(a); fmpz_init(b); fmpz_init(g);
    fmpz_init(k); fmpz_init(t); fmpz_init(x); fmpz_init(y);

    fmpz_divexact(k, c->m, c->q);

    /* a generator of the twists: a non-square (and non-cube for D = -3) */
    fmpz_sub_ui(t, n, 1);
    while (1)
    {
        fmpz_randm(g, state, n);
        if (fmpz_jacobi(g, n) != -1)
            continue;
        if (c->D == -3)
        {
            fmpz_divexact_ui(t, t, 3);
            fmpz_powm(b, g, t, n);
            fmpz_mul_ui(t, t, 3);
            if (fmpz_is_one(b))
                continue;
        }
        break;
    }

    if (c->D == -3)
    {
        /* y^2 = x^3 + g^i, i = 0..5 */
        fmpz_zero(a);
        fmpz_one(b);
        ntwists = 6;
    }
    else if (c->D == -4)
    {
        /* y^2 = x^3 + g^i x, i = 0..3 */
        fmpz_one(a);
        fmpz_zero(b);
        ntwists = 4;
    }
    else
    {
        if (!_ecpp_class_poly_root(j, c, state, ctx))
            goto cleanup;
        _ecpp_curve_from_j(a, b, j, ctx);
        ntwists = 2;
    }

    if (ntwists == 2 && flint_get_num_threads() >= 2)
    {
        /* both twists at once: costs one more scalar multiplication on
           average, halves the wall time */
        _twist_args args;
        fmpz xs[2], ys[2], tas[2], tbs[2];
        int res[2];
        slong w;

        for (w = 0; w < 2; w++)
        {
            fmpz_init(xs + w); fmpz_init(ys + w);
            fmpz_init(tas + w); fmpz_init(tbs + w);
        }
        args.a = a; args.b = b; args.g = g; args.q = c->q; args.k = k; args.ctx = ctx;
        args.x = xs; args.y = ys; args.ta = tas; args.tb = tbs; args.res = res;
        args.seed = n_randlimb(state);
        {
            PROF_START(t2);
            flint_parallel_do(_twist_worker, &args, 2, 2, FLINT_PARALLEL_STRIDED);
            PROF_ADD(prof_points, t2);
#ifdef ECPP_PROFILE
            prof_twists += 2;
#endif
        }
        for (w = 0; w < 2; w++)
        {
            if (res[w] == -1)
                result = -1;
            else if (res[w] == 1 && result != -1)
            {
                fmpz_set(s->n, n);
                s->D = c->D;
                fmpz_set(s->a, tas + w);
                fmpz_set(s->b, tbs + w);
                fmpz_set(s->m, c->m);
                fmpz_set(s->q, c->q);
                fmpz_set(s->x, xs + w);
                fmpz_set(s->y, ys + w);
                result = 1;
            }
        }
        for (w = 0; w < 2; w++)
        {
            fmpz_clear(xs + w); fmpz_clear(ys + w);
            fmpz_clear(tas + w); fmpz_clear(tbs + w);
        }
        goto cleanup;
    }

    for (twist = 0; twist < ntwists; twist++)
    {
        int r;

        if (twist > 0)
        {
            if (c->D == -3)
                fmpz_mod_mul(b, b, g, ctx);
            else if (c->D == -4)
                fmpz_mod_mul(a, a, g, ctx);
            else
            {
                /* (a, b) -> (a g^2, b g^3) */
                fmpz_mod_mul(t, g, g, ctx);
                fmpz_mod_mul(a, a, t, ctx);
                fmpz_mod_mul(t, t, g, ctx);
                fmpz_mod_mul(b, b, t, ctx);
            }
        }

        {
            PROF_START(t2);
            r = _ecpp_find_point(x, y, a, b, c->q, k, state, ctx);
            PROF_ADD(prof_points, t2);
#ifdef ECPP_PROFILE
            prof_twists++;
#endif
        }
        if (r == -1)
        {
            result = -1;
            goto cleanup;
        }
        if (r == 1)
        {
            fmpz_set(s->n, n);
            s->D = c->D;
            fmpz_set(s->a, a);
            fmpz_set(s->b, b);
            fmpz_set(s->m, c->m);
            fmpz_set(s->q, c->q);
            fmpz_set(s->x, x);
            fmpz_set(s->y, y);
            result = 1;
            goto cleanup;
        }
    }

cleanup:
    fmpz_clear(j); fmpz_clear(a); fmpz_clear(b); fmpz_clear(g);
    fmpz_clear(k); fmpz_clear(t); fmpz_clear(x); fmpz_clear(y);

    return result;
}

static void
_cand_init(ecpp_cand_struct * c)
{
    slong i;
    fmpz_init(c->m);
    fmpz_init(c->q);
    for (i = 0; i <= ECPP_DISC_MAXFAC; i++)
        fmpz_init(c->sqrts + i);
}

static void
_cand_clear(ecpp_cand_struct * c)
{
    slong i;
    fmpz_clear(c->m);
    fmpz_clear(c->q);
    for (i = 0; i <= ECPP_DISC_MAXFAC; i++)
        fmpz_clear(c->sqrts + i);
}

static void
_cand_set(ecpp_cand_struct * r, const ecpp_cand_struct * c)
{
    slong i;
    r->D = c->D; r->h = c->h; r->o = c->o; r->cost = c->cost;
    r->use_tower = c->use_tower; r->Dtower = c->Dtower; r->veven = c->veven;
    r->g = c->g;
    for (i = 0; i <= ECPP_DISC_MAXFAC; i++)
        r->pstar[i] = c->pstar[i];
    fmpz_set(r->m, c->m);
    fmpz_set(r->q, c->q);
    for (i = 0; i <= ECPP_DISC_MAXFAC; i++)
        fmpz_set(r->sqrts + i, c->sqrts + i);
}

static void
_step_init(ecpp_step_struct * st)
{
    fmpz_init(st->n); fmpz_init(st->a); fmpz_init(st->b); fmpz_init(st->m);
    fmpz_init(st->q); fmpz_init(st->x); fmpz_init(st->y);
    st->D = 0;
}

static void
_step_clear(ecpp_step_struct * st)
{
    fmpz_clear(st->n); fmpz_clear(st->a); fmpz_clear(st->b); fmpz_clear(st->m);
    fmpz_clear(st->q); fmpz_clear(st->x); fmpz_clear(st->y);
}

static void
_step_swap(ecpp_step_struct * a, ecpp_step_struct * b)
{
    slong t = a->D; a->D = b->D; b->D = t;
    fmpz_swap(a->n, b->n); fmpz_swap(a->a, b->a); fmpz_swap(a->b, b->b);
    fmpz_swap(a->m, b->m); fmpz_swap(a->q, b->q); fmpz_swap(a->x, b->x); fmpz_swap(a->y, b->y);
}

/* the pending realisation, as a thread pool task */
static void
_pending_task(void * varg)
{
    ecpp_ctx_struct * ctx = (ecpp_ctx_struct *) varg;
    ctx->pending.result = _ecpp_realise(&ctx->pending.step, ctx->pending.n,
                        &ctx->pending.cand, ctx->pending.state, ctx->pending.mctx);
}

/*
    Starts the realisation of (n, cand) as the pending one, on a pool
    thread if available, else immediately. The modulus context is
    duplicated so that the caller's may be cleared.
*/
static void
_pending_start(ecpp_ctx_struct * ctx, const fmpz_t n, const ecpp_cand_struct * cand,
                                                        slong cert_index)
{
    thread_pool_handle * handles = NULL;
    slong got = 0;

    ctx->pending.active = 1;
    ctx->pending.cert_index = cert_index;
    fmpz_set(ctx->pending.n, n);
    _cand_set(&ctx->pending.cand, cand);
    ctx->pending.mctx = flint_malloc(sizeof(fmpz_mod_ctx_struct));
    fmpz_mod_ctx_init(ctx->pending.mctx, n);
    flint_rand_set_seed(ctx->pending.state, n_randlimb(ctx->state), n_randlimb(ctx->state));

    /* one worker (the argument of flint_request_threads is a limit on the
       number of threads including the caller) */
    if (ctx->use_threads)
        got = flint_request_threads(&handles, 2);
    if (got == 1)
    {
        ctx->pending.threaded = 1;
        ctx->pending.handle = handles[0];
        flint_free(handles);
        thread_pool_wake(global_thread_pool, ctx->pending.handle, 0, _pending_task, ctx);
    }
    else
    {
        if (handles != NULL)
            flint_give_back_threads(handles, got);
        ctx->pending.threaded = 0;
        _pending_task(ctx);
    }
}

/*
    Joins the pending realisation, records its result for its certificate
    index and, on success, puts the step into the certificate. Returns the
    result (1, 0, -1), or 1 if nothing was pending.
*/
static int
_pending_join(ecpp_ctx_struct * ctx, ecpp_cert_t cert)
{
    slong idx;
    if (!ctx->pending.active)
        return 1;
    if (ctx->pending.threaded)
    {
        thread_pool_wait(global_thread_pool, ctx->pending.handle);
        thread_pool_give_back(global_thread_pool, ctx->pending.handle);
    }
    fmpz_mod_ctx_clear(ctx->pending.mctx);
    flint_free(ctx->pending.mctx);
    ctx->pending.active = 0;
    idx = ctx->pending.cert_index;
    if (idx >= ctx->results_alloc)
    {
        slong old = ctx->results_alloc;
        ctx->results_alloc = FLINT_MAX(2 * old, idx + 16);
        ctx->results = flint_realloc(ctx->results, ctx->results_alloc * sizeof(int));
        while (old < ctx->results_alloc)
            ctx->results[old++] = 0;
    }
    ctx->results[idx] = ctx->pending.result;
    if (ctx->pending.result == 1 && idx < cert->num)
        _step_swap(cert->steps + idx, &ctx->pending.step);
    return ctx->pending.result;
}

/*
    1: proved prime (steps appended), 0: composite, -1: gave up,
    -2: the realisation of the parent's step failed (its subtree,
    including this call, is abandoned).
*/
static int
_ecpp_prove_rec(ecpp_cert_t cert, const fmpz_t n, ecpp_ctx_struct * ctx, slong depth)
{
    ecpp_scan_struct scan;
    ecpp_cand_struct cand;
    fmpz_mod_ctx_t mctx;
    slong tries, my_index;
    int result = -1;

    /* fmpz_is_prime does not use ECPP for such small inputs */
    if (fmpz_bits(n) <= 64)
        return (fmpz_sgn(n) > 0 && fmpz_is_prime(n)) ? 1 : 0;

    if (depth > 4 * (slong) fmpz_bits(n))
        return -1;

    if (fmpz_is_even(n) || fmpz_divisible_si(n, 3))
        return 0;

    fmpz_mod_ctx_init(mctx, n);
    scan_init(&scan, n, mctx, ctx);
    _cand_init(&cand);

    for (tries = 0; tries < ctx->max_tries && result == -1; tries++)
    {
        int r, got, parent_ok;

        {
            PROF_START(t3);
            got = _scan_next_candidate(&scan, ctx, &cand);
            PROF_ADD(prof_scan, t3);
        }
        if (!got)
            break;

        /* before starting our realisation, join the parent's */
        parent_ok = _pending_join(ctx, cert);
        if (parent_ok != 1)
        {
            result = -2;
            break;
        }

        my_index = cert->num;
        ecpp_cert_push(cert);
        _pending_start(ctx, n, &cand, my_index);

#ifdef ECPP_PROFILE
        prof_nsteps++; prof_hsum += cand.o;
        prof_gain += fmpz_bits(n) - fmpz_bits(cand.q);
#endif
        if (ecpp_verbose)
            flint_printf("ecpp: %wd bits -> %wd bits, D = %wd (h = %wd, degree %wd)\n",
                fmpz_bits(n), fmpz_bits(cand.q), cand.D, cand.h, cand.o);
        /* not s->q: pushes may reallocate the steps */
        r = _ecpp_prove_rec(cert, cand.q, ctx, depth + 1);

        /* our realisation: joined by the child, or here */
        if (ctx->pending.active && ctx->pending.cert_index == my_index)
            _pending_join(ctx, cert);
        {
            int mine = ctx->results[my_index];
            if (mine != 1 && ecpp_verbose)
                flint_printf("ecpp: realisation failed (%d): D = %wd h = %wd\n", mine, cand.D, cand.h);
            if (mine == -1)
            {
                cert->num = my_index;
                result = 0;         /* composite */
                break;
            }
            if (mine != 1)
            {
                /* our step is invalid: the subtree below is discarded */
                cert->num = my_index;
                fmpz_set(scan.tried + scan.ntried++, cand.q);
                continue;
            }
        }
        if (r == 1)
        {
            result = 1;
            break;
        }
        /* composite q (impossible for a genuine candidate), dead end below,
           or abandoned by a failure above: next candidate */
        cert->num = my_index;
        fmpz_set(scan.tried + scan.ntried++, cand.q);
    }

    _cand_clear(&cand);
    scan_clear(&scan, ctx);
    fmpz_mod_ctx_clear(mctx);

    return result;
}

/*
    Returns 1 if n is proved prime (cert holds the certificate), 0 if n is
    proved composite, -1 if no proof could be found. n should be a probable
    prime; composites are usually detected quickly but not always cheaply.
*/
int
ecpp_prove(ecpp_cert_t cert, const fmpz_t n)
{
    ecpp_ctx_struct ctx;
    slong bits, i, pmax, hmax, tdivexp, Dmax, attempt;
    int result, bigpool;
    n_primes_t iter;

    cert->num = 0;

    if (fmpz_cmp_ui(n, 1) <= 0)
        return 0;
    if (fmpz_bits(n) <= 64)
        return fmpz_is_prime(n);                /* n > 1 here; no ECPP inside */
    if (!fmpz_is_probabprime(n))
        return 0;

    bits = fmpz_bits(n);
    _ecpp_tune(bits, &pmax, &hmax, &tdivexp);
    result = -1;

    /* if the pool is exhausted somewhere in the chain (rare, and only for
       small n), retry with a larger pool */
    /*
        If the pool is exhausted somewhere in the chain, retry with a
        larger one: more primes, and (the more important part) larger
        class numbers and odd parts admitted, at the price of more
        expensive realisations. This happens for the occasional n with
        few usable discriminants (e.g. n = 7 mod 8, which rules out every
        D = -4m and -8m with m = 3 mod 4) and no prime cofactor among them.
    */
    for (attempt = 0; result == -1 && pmax <= 100000; attempt++, pmax *= 2, hmax *= 2)
    {

    /* (the prime set is not reduced for the CM-style pool: the rounds
       take the primes in increasing order anyway, and a set that is too
       small leaves some n in the chain without candidates) */

    /* odd primes up to pmax */
    {
        slong alloc = 64;
        ulong p;
        ctx.primes = flint_malloc(alloc * sizeof(ulong));
        ctx.nprimes = 0;
        n_primes_init(iter);
        n_primes_next(iter);    /* skip 2 */
        while ((p = n_primes_next(iter)) <= (ulong) pmax)
        {
            if (ctx.nprimes == alloc)
            {
                alloc *= 2;
                ctx.primes = flint_realloc(ctx.primes, alloc * sizeof(ulong));
            }
            ctx.primes[ctx.nprimes++] = p;
        }
        n_primes_clear(iter);
    }

    /* the pool: smooth discriminants up to pmax^2, with a generous class
       number bound so that the search can persevere when the cheap part
       of the pool is exhausted */
    /*
        The pool: smooth discriminants up to Dmax. Only the odd part o of
        the class number matters for the realisation (the degree of the
        genus factor), so large class numbers are admitted as long as o is
        small (and h is not so large that the class polynomial itself
        becomes expensive); products of many small primes with |D| in the
        millions typically have h in the hundreds and o below 10. The
        table costs O(Dmax^{3/2}), hence the size dependent bound.
    */
    /*
        From 2500 bits on, the CM-style pool: all smooth discriminants
        up to Dmax with a cheap class field tower
        (measured at 4000 bits: 123 s instead of 137 s, square roots
        77 s -> 13 s; at 3000 bits parity, at 2000 bits slightly worse).
    */
    bigpool = (bits >= ECPP_BIGPOOL_BITS);
    Dmax = bigpool ? (WORD(1) << ECPP_DMAX_BITS(bits))
                   : FLINT_MIN(pmax * pmax, WORD(1) << ECPP_DMAX_BITS(bits));
    /*
        Not done yet: the class field tower with the Weber invariant makes
        discriminants with class numbers in the hundreds (h a product of
        small primes) realisable in about 0.1 s, as in Enge's CM library,
        which uses all such D up to 2^23 over a set of some 35 primes and
        gets ~45 bits per step from ~2500 Cornacchia attempts per step.
        Admitting them here (-DECPP_DISC_COST_MAX=60 and a larger Dmax)
        currently lowers the gain per step: the batch of orders is of
        fixed size, so more available orders do not become more candidates
        per probable prime test. Exploiting the larger pool needs the
        batch policy (size and trial division bound) retuned with it.
    */
    Dmax = FLINT_MAX(Dmax, 40000);
    ctx.ndiscs = ecpp_disc_table(&ctx.discs, ctx.primes, ctx.nprimes, Dmax,
                                        FLINT_MIN(4000, ECPP_POOL_HMAX << attempt),
                                        FLINT_MIN(256, ECPP_POOL_OMAX << attempt),
                                        (bigpool ? ECPP_DISC_COST_MAX : 0.0) * (attempt == 0 ? 1 : 3 << attempt));

    for (i = 0; i < 32; i++)
        ctx.have_primorial[i] = 0;
    ctx.max_tries = 16;
    ctx.tdivexp = tdivexp;
    flint_rand_init(ctx.state);
    ctx.pending.active = 0;
    ctx.pending.cert_index = -1;
    ctx.results_alloc = 8;
    ctx.results = flint_calloc(ctx.results_alloc, sizeof(int));
    _step_init(&ctx.pending.step);
    _cand_init(&ctx.pending.cand);
    fmpz_init(ctx.pending.n);
    flint_rand_init(ctx.pending.state);
    /* below a thousand bits or so the steps are too small for the thread
       pool's latency */
    ctx.use_threads = (flint_get_num_threads() >= 2 && bits >= 1000);

    if (ecpp_verbose)
        flint_printf("ecpp: %wd bits, %wd primes up to %wd, %wd discriminants up to %wd, trial division to 2^%wd\n",
                bits, ctx.nprimes, pmax, ctx.ndiscs, Dmax, tdivexp);

    result = _ecpp_prove_rec(cert, n, &ctx, 0);
    if (ctx.pending.active)
        _pending_join(&ctx, cert);
    if (result == -2)
        result = -1;

#ifdef ECPP_PROFILE
    flint_printf("ecpp profile: scan %.4fs classpoly %.4fs roots %.4fs points %.4fs | steps %wd avg degree %.1f avg gain %.1f twists tried %wd\n",
        prof_scan, prof_hilbert, prof_roots, prof_points, prof_nsteps,
        (double) prof_hsum / FLINT_MAX(prof_nsteps, 1), (double) prof_gain / FLINT_MAX(prof_nsteps, 1), prof_twists);
    flint_printf("  batches without a prime: %wd (%wd tests), batches tested to the end: %wd\n", prof_nobatch, prof_nobatch_tests, prof_fullbatch);
    {
        extern double ecpp_tower_prof[4];
        extern slong ecpp_tower_count, ecpp_tower_hsum, ecpp_tower_retries;
        flint_printf("  towers: %wd (avg h %.1f, %wd retries): class group %.4fs, values (arb) %.4fs, decomposition (arb) %.4fs, descent (mod n) %.4fs\n",
            ecpp_tower_count, (double) ecpp_tower_hsum / FLINT_MAX(1, ecpp_tower_count), ecpp_tower_retries,
            ecpp_tower_prof[0], ecpp_tower_prof[1], ecpp_tower_prof[2], ecpp_tower_prof[3]);
    }
    flint_printf("  scan: sqrt %.4fs (%wd) cornacchia %.4fs (%wd disc, %wd tried, %wd ok) split %.4fs (%wd orders, %wd batches) bpsw %.4fs (%wd)\n",
        prof_sqrt, prof_nsqrt, prof_corn, prof_ndisc, prof_ncorn, prof_ncornok, prof_split, prof_nsplit, prof_nbatch, prof_bpsw, prof_nbpsw);
#endif

    /* a proof must be a chain from n down: the steps are installed by
       whichever level joins their realisation, so check the structure
       (cheap; the verifier checks the mathematics) */
    if (result == 1)
    {
        slong k;
        int chain = (cert->num >= 1) && fmpz_equal(cert->steps[0].n, n);
        for (k = 1; k < cert->num && chain; k++)
            chain = fmpz_equal(cert->steps[k].n, cert->steps[k - 1].q);
        if (chain && fmpz_bits(cert->steps[cert->num - 1].q) > 64)
            chain = 0;
        if (!chain)
        {
            flint_printf("ecpp_prove: inconsistent certificate (internal error)\n");
            result = -1;
        }
    }
    if (result != 1)
        cert->num = 0;

    flint_rand_clear(ctx.state);
    _step_clear(&ctx.pending.step);
    _cand_clear(&ctx.pending.cand);
    flint_free(ctx.results);
    fmpz_clear(ctx.pending.n);
    flint_rand_clear(ctx.pending.state);
    for (i = 0; i < 32; i++)
        if (ctx.have_primorial[i])
        {
            slong k;
            fmpz_clear(ctx.primorial + i);
            if (ctx.use_threads)
                for (k = 0; k < ECPP_PRIMORIAL_CHUNKS; k++)
                    fmpz_clear(ctx.primorial_chunk[i] + k);
        }
    flint_free(ctx.primes);
    flint_free(ctx.discs);
    }

    return result;
}

int
ecpp_is_prime(const fmpz_t n)
{
    ecpp_cert_t cert;
    int r;

    ecpp_cert_init(cert);
    r = ecpp_prove(cert, n);
    ecpp_cert_clear(cert);

    return r;
}
