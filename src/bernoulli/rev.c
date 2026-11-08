/*
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    The reverse Bernoulli iterator.

    B_n for n = nmax, nmax - 2, ..., 0 is produced in batches of k
    consecutive indices n0, n0 - 2, ..., n0 - 2 (k - 1) (k = 16 for n0 >=
    4000, else 8) and yielded one by one from a buffer.

    A batch:

    1. The table pass: S_t = sum_(j odd >= 3) j^-(n0 - 2t) for t < k, from
       the stored powers j^-n0 in fixed point (frame 2^-wp), updated by
       multiplication by j^2.  The k steps are fused: each power is
       streamed once per batch.  The frame and the number of terms are
       lowered at each batch.

    2. The outputs: N_t = round(D_t P_t (1 + delta_t)), where D_t is the
       denominator, P_t = 2 n_t! / (2 pi)^(n_t) the prefactor (W limbs,
       floating) and delta_t = zeta(n_t) - 1 = S + (2^wp + S) / (2^n - 1)
       (the division by 2^n - 1 in linear time by suffix sums of n-bit
       chunks).  The rounding is certified by an error bound, with
       _bernoulli_fmpq_ui as a fallback.  Both products by P_t read the
       window [W - 1, 2 W) of a product of W-limb operands: delta is
       aligned to make this so, and with fft_small one plan serves the
       product precision for many batches (W is lowered in steps of 3%).

       Serially, the chain P_(t+1) = P_t (2 pi)^2 / (n (n - 1)), the
       transform of P_t shared by its two products and that of (2 pi)^2
       cached: four transforms per number.

       With threads, the jump P_t = P_0 (2 pi)^(2t) / (n0)_(2t) from the
       transform of P_0 and cached transforms of (2 pi)^(2t), which makes
       the numbers of a batch independent tasks (one division by a
       multi-limb integer more per number).  The threads take these tasks
       first (the expensive ones first), then chunks of powers of the
       table pass of the next batch, each thread summing into its own
       accumulators (two sets, by batch parity); the task of number t
       sums the accumulators of S_t.  The threads are requested when a
       batch starts, and the serial chain is used if none are available.

    The sieve (n >= SIEVE_CUTOFF): the table only holds the odd j prime
    to 3 and 5 (in the pass) and the {3, 5}-smooth j (summed apart into
    G, in limbs): 1 + sum_(j odd >= 3) j^-n = (1 + G)(1 + S) up to terms
    of distinct odd j > J, so S + G + G S with one short product per
    number; this saves about 7/15 of the table and of the pass.

    The table kernels: with AVX-512 IFMA, digits of 50 bits (below 2^51,
    redundant) in 64-bit lanes, valid for j < 2^24; otherwise full limbs
    with mpn_add_n and mpn_mul_1 (j^2 applied as two multiplications if it
    does not fit a limb), eight powers at a time by blocks of TILE limbs
    (each block stepped k times while in cache, the blocks of the k sums
    shared by the eight powers, the carries between blocks kept).

    Limits: n <= 2^30 (2^24 on 32-bit machines), so that bit counts fit a
    slong and the integers n (n - 1), j^2 and (n0)_(2t) fit a few limbs.
*/

#include <string.h>
#include <math.h>
#include "thread_pool.h"
#include "thread_support.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "arith.h"
#include "bernoulli.h"
#include "bernoulli/impl.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif
#if FLINT_USES_PTHREAD
# include <pthread.h>
#endif

#if FLINT_BITS == 64 && defined(__AVX512F__) && defined(__AVX512IFMA__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define REV_HAVE_IFMA 1
# include <immintrin.h>
#else
# define REV_HAVE_IFMA 0
#endif

#define REV_MODE_GENERIC 0
#define REV_MODE_IFMA 1

#define KMAX 16                 /* batch length */
#define GUARD 64                /* guard bits of the frame */
#define CHUNK 4096              /* powers between flushes of the digit sums */
#define ROOM (2 * KMAX + 16)    /* room for growth of the powers (digits) */
#define FFT_CUTOFF 150          /* fft_small products from W limbs */
#define THREAD_CUTOFF 400       /* threads from W limbs */
#define TABLE_GRAB 8            /* powers per grab in the threaded pass */
#define TILE 256                /* limbs per block of the tiled kernel */
#define SIEVE_CUTOFF 100000     /* the sieve by {3, 5} from this n */

/* ------------------------------------------------------------------------ */
/* helpers */

static double
rev_mag_2exp(const mag_t r, slong e)
{
    mag_t t;
    double d;
    mag_init(t);
    mag_mul_2exp_si(t, r, e);
    d = mag_get_d(t);
    mag_clear(t);
    return d;
}

/* x 2^e for any slong e */
static double
rev_ldexp(double x, slong e)
{
    e = FLINT_MAX(e, -100000);
    e = FLINT_MIN(e, 100000);
    return ldexp(x, (int) e);
}

static double
rev_log2(double x)
{
    return log(x) * 1.4426950408889634074;
}

/* the length of (x, n) without the high zero limbs, at least 1 */
static slong
rev_len(nn_srcptr x, slong n)
{
    while (n > 1 && x[n - 1] == 0)
        n--;
    return n;
}

static int
rev_batch_k(ulong n)
{
    return FLINT_MIN((n >= 4000) ? KMAX : 8, (n - BERNOULLI_REV_MIN) / 2 + 1);
}

static slong
rev_prod_limbs(slong wp)
{
    return (wp + 64 + FLINT_BITS - 1) / FLINT_BITS + 1;
}

/* x = m 2^e (m with exactly W limbs, top bit set) and the error in units
   of 2^e, from an arb ball */
static void
rev_arb_get_mpn(nn_ptr m, slong * e, double * err, const arb_t x, slong W)
{
    fmpz_t man, ex;
    slong b, sh;

    fmpz_init(man);
    fmpz_init(ex);
    arf_get_fmpz_2exp(man, ex, arb_midref(x));
    b = fmpz_bits(man);
    sh = FLINT_BITS * W - b;
    if (sh >= 0)
        fmpz_mul_2exp(man, man, sh);
    else
        fmpz_fdiv_q_2exp(man, man, -sh);
    *e = fmpz_get_si(ex) - sh;
    flint_mpn_zero(m, W);
    if (COEFF_IS_MPZ(*man))
    {
        mpz_ptr z = COEFF_TO_PTR(*man);
        flint_mpn_copyi(m, z->_mp_d, z->_mp_size);
    }
    else
        m[0] = *man;
    *err = rev_mag_2exp(arb_radref(x), -(*e)) * 1.000001 + (sh < 0 ? 1.0 : 0.0);
    fmpz_clear(man);
    fmpz_clear(ex);
}

/* a counter shared by the threads of a region */
typedef struct
{
    slong v;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
}
rev_counter_struct;

static void
rev_counter_init(rev_counter_struct * c)
{
    c->v = 0;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&c->mutex, NULL);
#endif
}

static void
rev_counter_clear(rev_counter_struct * c)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&c->mutex);
#else
    (void) c;
#endif
}

static slong
rev_counter_fetch_add(rev_counter_struct * c, slong d)
{
    slong r;
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&c->mutex);
#endif
    r = c->v;
    c->v += d;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&c->mutex);
#endif
    return r;
}

/* (z, zn) = floor(X Y / 2^wp) from the top limbs of X and Y (limbs above
   the frame, plus guard limbs: the result has at most tl limbs, computed
   from the lengths if tl = 0); returns the error in units, X and Y having
   errors ex and ey; z has room for xn + yn limbs */
static double
rev_mul_frame(nn_ptr z, slong * zn, nn_srcptr x, slong xn, double ex,
    nn_srcptr y, slong yn, double ey, slong wp, slong tl)
{
    slong a, b, sx, sy, sh, len;
    double e;

    if (tl == 0)
        tl = xn + yn - wp / FLINT_BITS + 1;
    tl = FLINT_MAX(tl, 1) + 3;
    a = FLINT_MIN(xn, tl);
    b = FLINT_MIN(yn, tl);
    sx = xn - a;
    sy = yn - b;
    sh = wp - FLINT_BITS * (sx + sy);
    if (sh < 0)
    {
        a = xn; b = yn; sx = sy = 0;
        sh = wp;
    }

    len = 0;
    if (a > 0 && b > 0 && sh / FLINT_BITS < a + b)
    {
        slong ls = sh / FLINT_BITS, bs = sh % FLINT_BITS;
        if (a >= b)
            flint_mpn_mul(z, x + sx, a, y + sy, b);
        else
            flint_mpn_mul(z, y + sy, b, x + sx, a);
        len = a + b - ls;
        if (bs != 0)
            mpn_rshift(z, z + ls, len, bs);
        else
            flint_mpn_copyi(z, z + ls, len);
        while (len > 0 && z[len - 1] == 0)
            len--;
    }
    *zn = len;

    /* the floor; the truncated limbs (X < 2^(B xn), the dropped parts
       below 2^(B sx) and 2^(B sy)); the errors of X and Y */
    e = 1.0;
    if (sy > 0)
        e += rev_ldexp(1.0, FLINT_BITS * (xn + sy) - wp);
    if (sx > 0)
        e += rev_ldexp(1.0, FLINT_BITS * (yn + sx) - wp);
    if (sx > 0 && sy > 0)
        e += rev_ldexp(1.0, FLINT_BITS * (sx + sy) - wp);
    e += rev_ldexp(ey, FLINT_BITS * xn - wp);
    e += rev_ldexp(ex, FLINT_BITS * yn - wp);
    e += rev_ldexp(ex * ey, -wp);
    return e * 1.000001;
}

/* whether j has no prime factors but those of the sieve (3, or 3 and 5) */
static int
rev_j_smooth(ulong j, int sieve)
{
    while (j % 3 == 0)
        j /= 3;
    if (sieve >= 2)
        while (j % 5 == 0)
            j /= 5;
    return j == 1;
}

/* whether j is prime to the sieve primes */
static int
rev_j_rough(ulong j, int sieve)
{
    return (sieve == 0) || (j % 3 != 0 && (sieve < 2 || j % 5 != 0));
}

/* ------------------------------------------------------------------------ */
/* table kernels: S_t += p c^t (t < k) and p *= c^k for c = j^2 */

/* limbs, S_t normalized (length slen), p with room for 2k more limbs; j^2
   is applied as two multiplications if it does not fit a limb (or if
   split is set) */
static slong
rev_kernel_generic(ulong * p, slong nd, ulong jj, int k, ulong ** S, slong slen,
    int split)
{
    int t;
    ulong cy, c = jj * jj;

    split = split || ((jj >> (FLINT_BITS / 2)) != 0);

    for (t = 0; t < k; t++)
    {
        if (nd == 0)
            break;
        cy = mpn_add_n(S[t], S[t], p, nd);
        if (cy)
            mpn_add_1(S[t] + nd, S[t] + nd, slen - nd, cy);
        if (!split)
        {
            cy = mpn_mul_1(p, p, nd, c);
            p[nd] = cy;
            nd += (cy != 0);
        }
        else
        {
            cy = mpn_mul_1(p, p, nd, jj);
            p[nd] = cy;
            nd += (cy != 0);
            cy = mpn_mul_1(p, p, nd, jj);
            p[nd] = cy;
            nd += (cy != 0);
        }
    }
    return nd;
}

#if REV_HAVE_IFMA

/* 64-bit limbs -> radix 2^db digits; returns the digit count */
static slong
rev_to_digits(ulong * d, nn_srcptr a, slong an, int db)
{
    slong nd = (an * 64 + db - 1) / db, i;
    ulong mask = (UWORD(1) << db) - 1;
    for (i = 0; i < nd; i++)
    {
        slong b = i * db, w = b / 64, s = b % 64;
        ulong v = a[w] >> s;
        if (s + db > 64 && w + 1 < an)
            v |= a[w + 1] << (64 - s);
        d[i] = v & mask;
    }
    while (nd > 0 && d[nd - 1] == 0)
        nd--;
    return nd;
}

/* (z, zn) += sum d_i 2^(db i), the digits below 2^63; x is scratch of zn
   limbs */
static void
rev_add_digits(nn_ptr z, slong zn, const ulong * d, slong nd, nn_ptr x, int db)
{
    slong i;
    ulong c = 0, mask = (UWORD(1) << db) - 1;
    flint_mpn_zero(x, zn);
    for (i = 0; i < nd || c != 0; i++)
    {
        ulong v = ((i < nd) ? d[i] : 0) + c;
        slong b = i * db, w = b / 64, s = b % 64;
        c = v >> db;
        v &= mask;
        x[w] |= v << s;
        if (s > 64 - db)
            x[w + 1] |= v >> (64 - s);
    }
    mpn_add_n(z, z, x, zn);
}

/* digits below 2^51, c < 2^48: with lo = 4 (x c mod 2^50) and
   hi = floor(4 x c / 2^52) = floor(x c / 2^50) < 2^49, the new digit
   lo / 4 + hi stays below 2^51 */
static inline __attribute__((always_inline)) slong
rev_kernel_ifma_k(ulong * p, slong nd, ulong c, const int k, ulong ** S)
{
    __m512i C4 = _mm512_set1_epi64(4 * c), Z = _mm512_setzero_si512();
    __m512i prev[KMAX];
    slong nb = (nd + k + 7) / 8, b, i;
    int t;

    _Pragma("GCC unroll 16")
    for (t = 0; t < k; t++)
        prev[t] = Z;
    for (i = nd; i < nb * 8; i++)
        p[i] = 0;

    for (b = 0; b < nb; b++)
    {
        __m512i x = _mm512_loadu_si512(p + 8 * b);
        __m512i s = _mm512_loadu_si512(S[0] + 8 * b);
        _mm512_storeu_si512(S[0] + 8 * b, _mm512_add_epi64(s, x));
        _Pragma("GCC unroll 16")
        for (t = 1; t <= k; t++)
        {
            __m512i lo = _mm512_madd52lo_epu64(Z, x, C4);
            __m512i hi = _mm512_madd52hi_epu64(Z, x, C4);
            lo = _mm512_srli_epi64(lo, 2);
            x = _mm512_add_epi64(lo, _mm512_alignr_epi64(hi, prev[t - 1], 7));
            prev[t - 1] = hi;
            if (t < k)
            {
                s = _mm512_loadu_si512(S[t] + 8 * b);
                _mm512_storeu_si512(S[t] + 8 * b, _mm512_add_epi64(s, x));
            }
        }
        _mm512_storeu_si512(p + 8 * b, x);
    }

    nd = nb * 8;
    while (nd > 0 && p[nd - 1] == 0)
        nd--;
    return nd;
}

static slong
rev_kernel_ifma(ulong * p, slong nd, ulong c, int k, ulong ** S)
{
    switch (k)
    {
        case 1: return rev_kernel_ifma_k(p, nd, c, 1, S);
        case 2: return rev_kernel_ifma_k(p, nd, c, 2, S);
        case 3: return rev_kernel_ifma_k(p, nd, c, 3, S);
        case 4: return rev_kernel_ifma_k(p, nd, c, 4, S);
        case 5: return rev_kernel_ifma_k(p, nd, c, 5, S);
        case 6: return rev_kernel_ifma_k(p, nd, c, 6, S);
        case 7: return rev_kernel_ifma_k(p, nd, c, 7, S);
        case 8: return rev_kernel_ifma_k(p, nd, c, 8, S);
        case 9: return rev_kernel_ifma_k(p, nd, c, 9, S);
        case 10: return rev_kernel_ifma_k(p, nd, c, 10, S);
        case 11: return rev_kernel_ifma_k(p, nd, c, 11, S);
        case 12: return rev_kernel_ifma_k(p, nd, c, 12, S);
        case 13: return rev_kernel_ifma_k(p, nd, c, 13, S);
        case 14: return rev_kernel_ifma_k(p, nd, c, 14, S);
        case 15: return rev_kernel_ifma_k(p, nd, c, 15, S);
        default: return rev_kernel_ifma_k(p, nd, c, 16, S);
    }
}

#endif

/* ------------------------------------------------------------------------ */
/* the state */

/* the parameters of a table pass */
typedef struct
{
    ulong n0;
    int k;
    slong wd;           /* frame: 2^-(db wd) */
    ulong J;            /* the largest j */
    slong sdlen;        /* digits of the sums */
    slong Sn;           /* limbs of the sums */
    int par;            /* the set of accumulators */
}
rev_batch_struct;

typedef struct
{
    /* table pass: digit sums (scratch) and the accumulators (2 sets) */
    ulong * Sd[KMAX];
    nn_ptr A[2][KMAX];
    double Aerr[2][KMAX];
    int touched[2];
    slong cnt;
    ulong car[2 * KMAX * TABLE_GRAB];  /* carries of the tiled kernel */
    /* outputs */
    nn_ptr tmp;
    slong tmpalloc;
#if FLINT_HAVE_FFT_SMALL
    slong gen;          /* the plan generation of the ops */
    fft_small_op_t X, Z;
#endif
}
rev_worker_struct;

typedef struct
{
    /* the output buffer */
    int nbuf, bpos;
    fmpz num[KMAX], den[KMAX];

    /* the power table */
    int mode;
    int split;          /* test: j^2 and n (n - 1) as two factors */
    int db;             /* bits per digit */
    slong wd;           /* the current frame: 2^-(db wd) */
    slong np;           /* powers j^-n for j = 3, 5, ..., 2 np + 1 */
    ulong ** buf;
    slong * off;
    slong * nd;
    slong * cap;
    double * err;       /* error bounds of the powers, frame units */
    slong Salloc, Sdalloc;

    /* the powers are of the odd j <= J listed in jv (increasing): with
       the sieve (by P = {3} or {3, 5}), only those prime to P (ri, the
       rough ones, nr of them, in the table pass) and the P-smooth ones
       (si, ns of them, kept in limbs and summed apart into G[par][t]) */
    int sieve;
    ulong J;
    ulong * jv;
    slong * ri, * si;
    slong nr, ns;
    nn_ptr G[2][KMAX];
    double Gerr[2][KMAX];

    /* the current batch and the next one (whose table pass is done when
       have_next is set) */
    rev_batch_struct cur, next;
    int have_next;
    int par;

    /* the sums of the current batch */
    nn_ptr S[KMAX];
    slong Sn[KMAX];
    double Serr[KMAX];

    /* prefactor 2 n! / (2 pi)^n = Pm 2^eP (W limbs), error Perr units 2^eP */
    nn_ptr Pm;
    slong W;
    slong eP;
    double Perr;

    /* (2 pi)^(2t) = ym[t] 2^ey[t] (Wy limbs), error yerr[t], t <= KMAX */
    nn_ptr ym[KMAX + 1];
    slong ey[KMAX + 1];
    double yerr[KMAX + 1];
    slong Wy;
    slong dmax;         /* limbs of the divisors (n0)_(2t) */

    /* workers */
    rev_worker_struct * wk;
    slong nwk;

#if FLINT_HAVE_FFT_SMALL
    int fft;
    slong gen;
    slong flo;
    ulong fslots;
    fft_small_plan_t FP;
    fft_small_op_struct Y[KMAX + 1];        /* transforms of y^t */
    int Yok[KMAX + 1];
    fft_small_op_t QH, QHn;                 /* of P_0, of the next P_0 */
    int QHok;                               /* QH holds the transform of Pm */
#endif
}
rev_struct;

typedef rev_struct rev_t[1];

static nn_ptr
rev_wk_tmp(rev_worker_struct * w, slong n)
{
    if (w->tmpalloc < n)
    {
        w->tmp = flint_realloc(w->tmp, n * sizeof(ulong));
        w->tmpalloc = n;
    }
    return w->tmp;
}

static void
rev_fit_workers(rev_t B, slong nw)
{
    slong i;
    if (nw <= B->nwk)
        return;
    B->wk = flint_realloc(B->wk, nw * sizeof(rev_worker_struct));
    for (i = B->nwk; i < nw; i++)
    {
        memset(B->wk + i, 0, sizeof(rev_worker_struct));
#if FLINT_HAVE_FFT_SMALL
        B->wk[i].gen = -1;
#endif
    }
    B->nwk = nw;
}

/* the table buffers of a worker, allocated when it first takes powers
   (the generic kernel sums into the accumulators directly) */
static void
rev_wk_fit_table(rev_t B, rev_worker_struct * w)
{
    int t, s;
    if (w->A[0][0] != NULL)
        return;
    for (t = 0; t < KMAX; t++)
    {
        if (B->mode != REV_MODE_GENERIC)
            w->Sd[t] = flint_aligned_alloc(64, ((B->Sdalloc + 8) * sizeof(ulong) + 63) / 64 * 64);
        for (s = 0; s < 2; s++)
            w->A[s][t] = flint_malloc(B->Salloc * sizeof(ulong));
    }
}

static void
rev_clear_workers(rev_t B)
{
    slong i;
    int t, s;
    for (i = 0; i < B->nwk; i++)
    {
        rev_worker_struct * w = B->wk + i;
        if (w->A[0][0] != NULL)
        {
            for (t = 0; t < KMAX; t++)
            {
                if (w->Sd[t] != NULL)
                    flint_aligned_free(w->Sd[t]);
                for (s = 0; s < 2; s++)
                    flint_free(w->A[s][t]);
            }
        }
        flint_free(w->tmp);
#if FLINT_HAVE_FFT_SMALL
        if (w->gen >= 0)
        {
            fft_small_op_clear(w->X);
            fft_small_op_clear(w->Z);
        }
#endif
    }
    flint_free(B->wk);
    B->wk = NULL;
    B->nwk = 0;
}

/* ------------------------------------------------------------------------ */
/* the fft plan for the products by Pm: window [flo, 2W) */

#if FLINT_HAVE_FFT_SMALL

static void
rev_fft_clear(rev_t B)
{
    int t;
    if (B->fft)
    {
        for (t = 1; t <= KMAX; t++)
            if (B->Yok[t])
                fft_small_op_clear(B->Y + t);
        fft_small_op_clear(B->QH);
        fft_small_op_clear(B->QHn);
        fft_small_plan_clear(B->FP);
        B->fft = 0;
        for (t = 1; t <= KMAX; t++)
            B->Yok[t] = 0;
        B->QHok = 0;
    }
}

static void
rev_fft_setup(rev_t B)
{
    slong W = B->W;
    rev_fft_clear(B);
    B->gen++;
    if (W < FFT_CUTOFF)
        return;
    B->flo = W - 2 - B->dmax;
    if (!fft_small_plan_init_mpn_window(B->FP, get_default_mpn_ctx(), W, W,
            1, 1, B->flo, 2 * W))
        return;
    B->fslots = fft_small_plan_mpn_nslots(B->FP, W, W);
    fft_small_op_init(B->QH, B->FP);
    fft_small_op_init(B->QHn, B->FP);
    B->QHok = 0;
    B->fft = 1;
}

/* the transform of y^t (computed by the only task that needs it) */
static void
rev_fft_fit_Y(rev_t B, int t)
{
    if (!B->Yok[t])
    {
        fft_small_op_init(B->Y + t, B->FP);
        fft_small_fft_mpn(B->Y + t, B->ym[t] + (B->Wy - B->W), B->W, B->FP);
        B->Yok[t] = 1;
    }
}

static void
rev_wk_fit_ops(rev_t B, rev_worker_struct * w)
{
    if (w->gen != B->gen)
    {
        if (w->gen >= 0)
        {
            fft_small_op_clear(w->X);
            fft_small_op_clear(w->Z);
        }
        fft_small_op_init(w->X, B->FP);
        fft_small_op_init(w->Z, B->FP);
        w->gen = B->gen;
    }
}

/* the window [flo, 2W) of the product in X, into z (2W - flo limbs) */
static void
rev_fft_export(nn_ptr z, rev_t B, fft_small_op_t X)
{
    fft_small_ifft_trunc(X, B->fslots, 0, NULL, B->FP);
    fft_small_export_mpn_signed_window(z, 2 * B->W - B->flo, X, B->fslots, B->flo, B->FP);
}

/* the same for a product of Pm by an operand of bn limbs: only the slots
   of the product length are inverse transformed and exported (the window
   [flo, W + bn + 1), zero padded to 2W - flo limbs) */
static void
rev_fft_export_short(nn_ptr z, rev_t B, fft_small_op_t X, slong bn)
{
    slong W = B->W, hi = FLINT_MIN(2 * W, W + bn + 1);
    ulong ns = FLINT_MIN(B->fslots, fft_small_plan_mpn_nslots(B->FP, W, bn) + 1);
    fft_small_ifft_trunc(X, ns, 0, NULL, B->FP);
    fft_small_export_mpn_signed_window(z, hi - B->flo, X, ns, B->flo, B->FP);
    flint_mpn_zero(z + (hi - B->flo), 2 * W - hi);
}

#endif

/* ------------------------------------------------------------------------ */
/* the table pass */

static slong
rev_frame_digits(rev_t B, ulong n)
{
    slong wp = bernoulli_global_prec(n) + GUARD;
    return (wp + B->db - 1) / B->db;
}

/* sum_(j odd > J) j^-m in units of 2^-wp: at most
   (J + 2)^-m + J^(1-m) / (2 (m - 1)) */
static double
rev_tail(ulong J, ulong m, slong wp)
{
    double a = wp - m * rev_log2(J + 2.0);
    double b = wp + (1.0 - m) * rev_log2(J) - rev_log2(2.0 * (m - 1));
    return (rev_ldexp(1.0, (slong) ceil(FLINT_MIN(a, 2000.0)))
          + rev_ldexp(1.0, (slong) ceil(FLINT_MIN(b, 2000.0)))) * 1.0001;
}

/* prepare the table pass of the batch at n0: lower the frame, drop the
   powers that are no longer needed */
static void
rev_prepare(rev_t B, rev_batch_struct * T, ulong n0)
{
    slong wd_new, d, i, np_new, wp, sdlen;

    T->n0 = n0;
    T->k = rev_batch_k(n0);
    T->par = B->par;
    B->par ^= 1;

    wd_new = rev_frame_digits(B, n0);
    if (wd_new < B->wd)
    {
        d = B->wd - wd_new;
        for (i = 0; i < B->np; i++)
        {
            if (B->mode != REV_MODE_GENERIC && B->sieve && !rev_j_rough(B->jv[i], B->sieve))
            {
                /* a smooth power in limbs: drop db d bits */
                slong bits = B->db * d, ls = bits / FLINT_BITS, bs = bits % FLINT_BITS;
                ls = FLINT_MIN(ls, B->nd[i]);
                B->off[i] += ls;
                B->nd[i] -= ls;
                if (bs != 0 && B->nd[i] > 0)
                {
                    mpn_rshift(B->buf[i] + B->off[i], B->buf[i] + B->off[i], B->nd[i], bs);
                    while (B->nd[i] > 0 && B->buf[i][B->off[i] + B->nd[i] - 1] == 0)
                        B->nd[i]--;
                }
            }
            else
            {
                slong dd = FLINT_MIN(d, B->nd[i]);
                B->off[i] += dd;
                B->nd[i] -= dd;
            }
            /* the dropped digits (each below 2^(db + 1)) are worth less
               than two units */
            B->err[i] = rev_ldexp(B->err[i], -B->db * d) * 1.000001 + 2.0;
        }
        B->wd = wd_new;
    }
    wp = B->wd * B->db;
    T->wd = B->wd;

    B->J = FLINT_MIN(B->J, bernoulli_zeta_terms(n0 - 2 * (T->k - 1), wp));
    np_new = B->np;
    while (np_new > 0 && B->jv[np_new - 1] > B->J)
        np_new--;
    for (i = np_new; i < B->np; i++)
        flint_aligned_free(B->buf[i]);
    B->np = np_new;
    while (B->nr > 0 && B->ri[B->nr - 1] >= np_new)
        B->nr--;
    while (B->ns > 0 && B->si[B->ns - 1] >= np_new)
        B->ns--;
    T->J = B->J;

    /* the digit sums of the rough powers; the limb sums of all, the
       smooth ones (in limbs) included */
    sdlen = 0;
    for (i = 0; i < B->nr; i++)
        sdlen = FLINT_MAX(sdlen, B->nd[B->ri[i]]);
    sdlen += ROOM;
    FLINT_ASSERT(sdlen <= B->Sdalloc);
    T->sdlen = sdlen;
    T->Sn = (sdlen * B->db + FLINT_BITS - 1) / FLINT_BITS + 2;
    for (i = 0; i < B->ns; i++)
        T->Sn = FLINT_MAX(T->Sn, B->nd[B->si[i]] + ROOM + 2);
    FLINT_ASSERT(T->Sn <= B->Salloc);
}

/* the generic kernel on the powers idx[0], ..., idx[cnt - 1] at once, by
   blocks of TILE limbs: each block of a power is stepped k times while
   in cache and the blocks of the k sums are reused by all the powers,
   instead of streaming each power and the k sums k times (which is
   bound by the memory bandwidth for powers beyond the L2 cache); the
   carries of the products between blocks are kept in car (2 k per power) */
static void
rev_kernel_generic_tiled(rev_struct * B, const slong * idx, slong cnt, int k,
    ulong ** S, slong slen, int split, ulong * car)
{
    slong top = 0, e, o, tops[TABLE_GRAB];
    ulong c[TABLE_GRAB];
    int two[TABLE_GRAB], t;

    for (e = 0; e < cnt; e++)
    {
        slong i = idx[e];
        ulong jj = B->jv[i];
        nn_ptr p = B->buf[i] + B->off[i];

        two[e] = split || ((jj >> (FLINT_BITS / 2)) != 0);
        c[e] = two[e] ? jj : jj * jj;
        tops[e] = (B->nd[i] == 0) ? 0 : B->nd[i] + (two[e] ? 2 : 1) * k;
        if (tops[e] > B->nd[i])
            flint_mpn_zero(p + B->nd[i], tops[e] - B->nd[i]);
        top = FLINT_MAX(top, tops[e]);
        for (t = 0; t < 2 * k; t++)
            car[2 * KMAX * e + t] = 0;
    }

    for (o = 0; o < top; o += TILE)
    {
        for (e = 0; e < cnt; e++)
        {
            slong len;
            nn_ptr p;
            ulong * cr = car + 2 * KMAX * e;

            if (o >= tops[e])
                continue;
            len = FLINT_MIN(TILE, tops[e] - o);
            p = B->buf[idx[e]] + B->off[idx[e]] + o;

            for (t = 0; t < k; t++)
            {
                ulong cy = mpn_add_n(S[t] + o, S[t] + o, p, len), hi;
                if (cy)
                {
                    FLINT_ASSERT(o + len < slen);
                    mpn_add_1(S[t] + o + len, S[t] + o + len, slen - o - len, cy);
                }
                hi = mpn_mul_1(p, p, len, c[e]);
                hi += mpn_add_1(p, p, len, cr[2 * t]);
                cr[2 * t] = hi;
                if (two[e])
                {
                    hi = mpn_mul_1(p, p, len, c[e]);
                    hi += mpn_add_1(p, p, len, cr[2 * t + 1]);
                    cr[2 * t + 1] = hi;
                }
            }
        }
    }

    for (e = 0; e < cnt; e++)
    {
        slong i = idx[e], nd = tops[e];
        nn_srcptr p = B->buf[i] + B->off[i];
        while (nd > 0 && p[nd - 1] == 0)
            nd--;
        B->nd[i] = nd;
    }
}

static void
rev_table_prep(rev_struct * B, slong i, const rev_batch_struct * T, double * Aerr)
{
    ulong jj = B->jv[i];
    double ce = 1.0, c = (double) jj * (double) jj;
    slong need;
    int t, k = T->k;

    /* a new buffer when the power outgrows its buffer or has shrunk
       well below it (the frame being lowered faster than the power
       grows, the memory is given back as the iteration proceeds) */
    need = B->nd[i] + ROOM;
    if (B->off[i] + need > B->cap[i] || B->cap[i] > need + need / 8 + 128)
    {
        if (need > B->cap[i] || B->cap[i] > need + need / 8 + 128)
        {
            slong ncap = need + 64;
            ulong * nb = flint_aligned_alloc(64, ((ncap * sizeof(ulong) + 63) / 64) * 64);
            flint_mpn_copyi(nb, B->buf[i] + B->off[i], B->nd[i]);
            flint_aligned_free(B->buf[i]);
            B->buf[i] = nb;
            B->cap[i] = ncap;
        }
        else
            memmove(B->buf[i], B->buf[i] + B->off[i], B->nd[i] * sizeof(ulong));
        B->off[i] = 0;
    }

    for (t = 0; t < k; t++)
    {
        Aerr[t] += B->err[i] * ce;
        ce *= c;
    }
    B->err[i] *= ce * 1.000001;
}

/* one power: grow its buffer if needed, accumulate its error terms into
   Aerr, run the kernel: into the limb sums S of length slen (a smooth
   power, S != NULL), or into the digit sums of the worker (IFMA; the
   rough powers of the generic mode go to rev_kernel_generic_tiled) */
static void
rev_table_one(rev_struct * B, rev_worker_struct * w, slong i, const rev_batch_struct * T,
    nn_ptr * S, slong slen, double * Aerr)
{
    ulong jj = B->jv[i];
    int k = T->k;

    rev_table_prep(B, i, T, Aerr);

    if (S != NULL)
    {
        B->nd[i] = rev_kernel_generic(B->buf[i] + B->off[i], B->nd[i], jj, k,
            S, slen, B->split);
        return;
    }
#if REV_HAVE_IFMA
    FLINT_ASSERT(B->mode == REV_MODE_IFMA);
    B->nd[i] = rev_kernel_ifma(B->buf[i] + B->off[i], B->nd[i], jj * jj, k, w->Sd);
#else
    (void) w;
    FLINT_ASSERT(0);
#endif
}

/* the smooth powers of the pass T, into G[par] */
static void
rev_smooth_pass(rev_struct * B, const rev_batch_struct * T)
{
    slong i;
    int t;

    for (t = 0; t < T->k; t++)
    {
        flint_mpn_zero(B->G[T->par][t], T->Sn);
        B->Gerr[T->par][t] = 0.0;
    }
    for (i = 0; i < B->ns; i++)
        rev_table_one(B, NULL, B->si[i], T, B->G[T->par], T->Sn, B->Gerr[T->par]);
}

/* the digit sums of a worker into its accumulators */
static void
rev_wk_flush(rev_struct * B, rev_worker_struct * w, const rev_batch_struct * T)
{
    if (B->mode == REV_MODE_GENERIC)
        return;
#if REV_HAVE_IFMA
    {
        nn_ptr x = rev_wk_tmp(w, T->Sn + 8);
        int t;
        for (t = 0; t < T->k; t++)
        {
            rev_add_digits(w->A[T->par][t], T->Sn, w->Sd[t], T->sdlen, x, B->db);
            flint_mpn_zero(w->Sd[t], T->sdlen);
        }
    }
#endif
    w->cnt = 0;
}

/* a worker starts taking powers for the pass T */
static void
rev_wk_table_start(rev_struct * B, rev_worker_struct * w, const rev_batch_struct * T)
{
    int t;
    rev_wk_fit_table(B, w);
    for (t = 0; t < T->k; t++)
    {
        /* the generic sums (in A) are limbs of length sdlen < Sn */
        if (B->mode != REV_MODE_GENERIC)
            flint_mpn_zero(w->Sd[t], FLINT_MAX(T->sdlen, T->Sn));
        flint_mpn_zero(w->A[T->par][t], T->Sn);
        w->Aerr[T->par][t] = 0.0;
    }
    w->cnt = 0;
    w->touched[T->par] = 1;
}

/* S_t of a batch: the sum of the accumulators of all workers */
static void
rev_reduce(rev_t B, const rev_batch_struct * T, int t)
{
    slong i, Sn = T->Sn;
    int first = 1;
    double e = 0.0;

    for (i = 0; i < B->nwk; i++)
    {
        rev_worker_struct * w = B->wk + i;
        if (!w->touched[T->par])
            continue;
        if (first)
            flint_mpn_copyi(B->S[t], w->A[T->par][t], Sn);
        else
            mpn_add_n(B->S[t], B->S[t], w->A[T->par][t], Sn);
        e += w->Aerr[T->par][t];
        first = 0;
    }
    if (first)
        flint_mpn_zero(B->S[t], Sn);
    B->Serr[t] = e * 1.000001;
    B->Sn[t] = rev_len(B->S[t], Sn);

    /* with the sieve: sum_(j odd) j^-n = (sum_(smooth j) j^-n)
       (sum_(rough j) j^-n), that is S + G + G S */
    if (B->sieve)
    {
        nn_srcptr G = B->G[T->par][t];
        slong gn = rev_len(G, Sn), wp = T->wd * B->db, zn;
        nn_ptr z = flint_malloc(2 * (Sn + 2) * sizeof(ulong));
        double ez;

        ez = rev_mul_frame(z, &zn, G, gn, B->Gerr[T->par][t], B->S[t], B->Sn[t],
            B->Serr[t], wp, 0);
        mpn_add_n(B->S[t], B->S[t], G, Sn);
        if (zn > 0)
            mpn_add(B->S[t], B->S[t], Sn, z, zn);
        B->Serr[t] += B->Gerr[T->par][t] * 1.000001 + ez;
        B->Sn[t] = rev_len(B->S[t], Sn);
        flint_free(z);
    }

    /* the omitted terms: with the sieve, those of the odd j whose smooth
       or rough part exceeds J (the products of a smooth and a rough
       j <= J are exact terms, whatever their size), so in any case odd
       j > J only, each at most once; the sums being one-sided, the true
       value lies within the tail above the computed one */
    B->Serr[t] += rev_tail(T->J, T->n0 - 2 * t, T->wd * B->db);
    while (B->Sn[t] > 0 && B->S[t][B->Sn[t] - 1] == 0)
        B->Sn[t]--;
}

/* ------------------------------------------------------------------------ */
/* the outputs */

/* a prefactor P = Pm 2^eP, W limbs, error Perr units */
typedef struct
{
    nn_ptr m;
    slong e;
    double err;
}
rev_pref_struct;

/* (z, zn) += floor(X / (2^n - 1)) approximately, X = (x, zn): with
   n-bit chunks x_k of X, the quotient is sum_m s_(m+1) 2^(m n) plus the
   fraction s_0 / (2^n - 1), s_m = sum_(k >= m) x_k; returns a bound for
   the error (the fraction is below the number of chunks + 1); the sum
   must fit in zn limbs; U has room for 2 (n / B + 3) + 2 limbs */
static double
rev_add_div_2expm1(nn_ptr z, nn_srcptr x, slong zn, ulong n, nn_ptr U)
{
    slong cl = n / FLINT_BITS + 3, nch, m;
    nn_ptr s = U, r = U + cl;
    slong xbits = zn * FLINT_BITS;

    while (xbits > 0 && x[(xbits - 1) / FLINT_BITS] == 0)
        xbits -= FLINT_BITS;
    nch = xbits / n + 1;
    flint_mpn_zero(s, cl);

    for (m = nch - 1; m >= 0; m--)
    {
        ulong b = (m + 1) * n, bl = b / FLINT_BITS, bb = b % FLINT_BITS;
        ulong o = m * n, ol = o / FLINT_BITS, ob = o % FLINT_BITS;
        slong cnt, i, an;

        /* s += chunk m + 1 */
        if ((slong) bl < zn)
        {
            cnt = FLINT_MIN(cl - 1, zn - (slong) bl);
            if (bb)
                mpn_rshift(r, x + bl, cnt, bb);
            else
                flint_mpn_copyi(r, x + bl, cnt);
            for (i = cnt; i < cl; i++)
                r[i] = 0;
            if (n % FLINT_BITS)
                r[n / FLINT_BITS] &= (UWORD(1) << (n % FLINT_BITS)) - 1;
            else
                r[n / FLINT_BITS] = 0;
            for (i = n / FLINT_BITS + 1; i < cl; i++)
                r[i] = 0;
            mpn_add_n(s, s, r, cl);
        }

        /* z += s 2^(m n) */
        if (ob)
            r[cl] = mpn_lshift(r, s, cl, ob);
        else
        {
            flint_mpn_copyi(r, s, cl);
            r[cl] = 0;
        }
        an = cl + 1;
        while (an > 0 && r[an - 1] == 0)
            an--;
        if (an == 0)
            continue;
        FLINT_ASSERT((slong) ol + an <= zn);
        mpn_add(z + ol, z + ol, zn - ol, r, an);
    }

    return nch + 2.0;
}

/* delta' = (zeta(n) - 1) 2^wp B^e with e = W - 1 - floor(wp / B), in
   (dl, W) (W + 2 limbs of room); returns the error in units of 2^-wp;
   scratch U has 2 W + 2 (n / B) + 16 limbs */
static double
rev_delta(nn_ptr dl, nn_srcptr S, slong Sn, double eS, ulong n, slong wp,
    slong W, nn_ptr U)
{
    slong s1 = wp / FLINT_BITS, e = W - 1 - s1, ln = s1 + 2;
    nn_ptr L = U, V = U + ln;
    double ed;

    FLINT_ASSERT(e >= 0 && Sn <= ln);
    flint_mpn_zero(L, ln);
    flint_mpn_copyi(L, S, Sn);
    L[s1] += UWORD(1) << (wp % FLINT_BITS);
    flint_mpn_zero(dl, W + 2);
    flint_mpn_copyi(dl + e, S, Sn);
    ed = eS + rev_add_div_2expm1(dl + e, L, ln, n, V);
    ed = ed * (1.0 + rev_ldexp(1.0, 1 - (slong) n)) + 2.0;
    FLINT_ASSERT(dl[W] == 0);
    return ed;
}

/* N = round(D P (1 + delta)) from T = floor(Pm delta / 2^wp), whose
   error is eTabs + eTrel 2^sigma units of 2^eP; the fraction of the
   rounded value must stay away from 1/2 by more than the error, else B_n
   is computed by _bernoulli_fmpq_ui */
static void
rev_round(fmpz_t numer, fmpz_t denom, ulong n, const rev_pref_struct * P,
    slong W, nn_srcptr Tm, slong Tn, double eTabs, double eTrel, slong sigma,
    const fmpz_t D, nn_ptr scr)
{
    slong ptn, Vn, sh;
    nn_ptr PT = scr, V = scr + W + Tn + 2;
    int ok;

    ptn = FLINT_MAX(W, Tn) + 1;
    flint_mpn_zero(PT, ptn);
    flint_mpn_copyi(PT, P->m, W);
    if (Tn > 0)
        mpn_add(PT, PT, ptn, Tm, Tn);
    while (ptn > 0 && PT[ptn - 1] == 0)
        ptn--;
    if (!COEFF_IS_MPZ(*D))
    {
        V[ptn] = mpn_mul_1(V, PT, ptn, *D);
        Vn = ptn + 1;
    }
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*D);
        flint_mpn_mul(V, PT, ptn, z->_mp_d, z->_mp_size);
        Vn = ptn + z->_mp_size;
    }
    V[Vn] = 0;
    Vn++;

    /* round at bit sh, checking the top limb of the fraction */
    sh = -P->e;
    FLINT_ASSERT(sh > FLINT_BITS);
    {
        slong sl = sh / FLINT_BITS, sb = sh % FLINT_BITS;
        slong b = sh - FLINT_BITS, bl = b / FLINT_BITS, bb = b % FLINT_BITS;
        ulong top;
        double q, lim = ldexp(1.0, FLINT_BITS);

        mpn_add_1(V + (sh - 1) / FLINT_BITS, V + (sh - 1) / FLINT_BITS,
            Vn - (sh - 1) / FLINT_BITS, UWORD(1) << ((sh - 1) % FLINT_BITS));
        top = V[bl] >> bb;
        if (bb)
            top |= V[bl + 1] << (FLINT_BITS - bb);
        /* the error in units of 2^(sh - B), plus a margin for the
           conversion of top to a double */
        q = fmpz_get_d(D) * 1.000001 * (rev_ldexp(P->err * 1.000001 + eTabs + 1.0, -b)
                + rev_ldexp(eTrel, sigma - b)) + ((FLINT_BITS == 64) ? 8192.0 : 1.0);
        ok = q < (double) top && (double) top + q < lim;

        if (ok)
        {
            nn_ptr Q = V + sl;
            slong qn = Vn - sl;
            if (sb)
                mpn_rshift(Q, Q, qn, sb);
            while (qn > 0 && Q[qn - 1] == 0)
                qn--;
            fmpz_set_ui_array(numer, Q, qn);
            if (n % 4 == 0)
                fmpz_neg(numer, numer);
            fmpz_set(denom, D);
        }
        else
        {
            _bernoulli_fmpq_ui(numer, denom, n);
        }
    }
}

/* scratch limbs for one output */
static slong
rev_out_scratch(slong W, ulong n, slong dmax)
{
    return 18 * W + 6 * dmax + 2 * (slong) (n / FLINT_BITS) + 256;
}

/* T from the window of Pm delta' (limbs [W-1, 2W) shifted by s0 bits) */
static slong
rev_T_from_window(nn_ptr Tm, nn_srcptr z, slong W, slong s0)
{
    slong tn = W + 1;
    if (s0)
        mpn_rshift(Tm, z, tn, s0);
    else
        flint_mpn_copyi(Tm, z, tn);
    while (tn > 0 && Tm[tn - 1] == 0)
        tn--;
    return tn;
}

/* Perr of a new prefactor from relative errors rel (units of B^-(W-1)) */
static double
rev_new_err(nn_srcptr m, slong W, double rel, double abs)
{
    return ((double) m[W - 1] + 1.0) * rel * 1.00001 + abs;
}

/* (x, xn) = floor(x / (n (n - 1))) */
static void
rev_div_nn1(nn_ptr x, slong xn, ulong n, int split)
{
    ulong hi, lo;
    umul_ppmm(hi, lo, n, n - 1);
    if (hi == 0 && !split)
        mpn_divrem_1(x, 0, x, xn, lo);
    else
    {
        mpn_divrem_1(x, 0, x, xn, n);
        mpn_divrem_1(x, 0, x, xn, n - 1);
    }
}

/* --- the serial chain: outputs t = 0, ..., k-1 and the next prefactor --- */

static void
rev_outputs_chain(rev_t B, rev_worker_struct * w)
{
    const rev_batch_struct * T = &B->cur;
    int k = T->k;
    ulong n0 = T->n0;
    slong W = B->W, wp = T->wd * B->db, s0 = wp % FLINT_BITS, t;
    slong Wy = B->Wy;
    nn_srcptr y = B->ym[1] + (Wy - W);
    slong ey = B->ey[1] + FLINT_BITS * (Wy - W);
    double yerr = B->yerr[1] + 1.0;

    for (t = 0; t < k; t++)
        rev_reduce(B, T, t);

    for (t = 0; t < k; t++)
    {
        ulong n = n0 - 2 * t;
        nn_ptr scr = rev_wk_tmp(w, rev_out_scratch(W, n, B->dmax));
        nn_ptr dl = scr, H = dl + W + 2, Tm = H + 2 * W + 2, z = Tm + W + 2;
        nn_ptr U = z + 2 * W + 2;
        slong Tn;
        double ed, rel;
        rev_pref_struct P;

        ed = rev_delta(dl, B->S[t], B->Sn[t], B->Serr[t], n, wp, W, U);

#if FLINT_HAVE_FFT_SMALL
        if (B->fft)
        {
            slong off = W - 1 - B->flo, dn;
            rev_wk_fit_ops(B, w);
            rev_fft_fit_Y(B, 1);
            fft_small_fft_mpn(w->X, B->Pm, W, B->FP);
            dn = rev_len(dl, W);
            fft_small_fft_mpn(w->Z, dl, dn, B->FP);
            fft_small_op_mul(w->Z, w->X, w->Z, B->FP);
            fft_small_op_mul(w->X, w->X, B->Y + 1, B->FP);
            rev_fft_export(z, B, w->X);
            flint_mpn_copyi(H, z + off, W + 1);
            rev_fft_export_short(z, B, w->Z, dn);
            Tn = rev_T_from_window(Tm, z + off, W, s0);
        }
        else
#endif
        {
            /* within 2W units of the windows (flint_mpn_mulhigh_n) */
            H[0] = flint_mpn_mulhigh_n(H + 1, B->Pm, y, W);
            z[0] = flint_mpn_mulhigh_n(z + 1, B->Pm, dl, W);
            Tn = rev_T_from_window(Tm, z, W, s0);
        }

        P.m = B->Pm;
        P.e = B->eP;
        P.err = B->Perr;
        rev_round(B->num + t, B->den + t, n, &P, W, Tm, Tn, 2.0 * W + 4.0,
            ((double) B->Pm[W - 1] + 1.0) * ed, FLINT_BITS * (W - 1) - wp,
            B->den + t, U);

        /* P <- P y / (n (n - 1)), normalized to a top bit set (as the
           initial P, so that the top limb stays nonzero: P y has its top
           bit in one of its top two limbs and n (n - 1) < 2^(2 B - 2));
           the bits shifted in below H[0] are zeros, and the error in the
           units of the new P grows by the shift */
        rev_div_nn1(H, W + 1, n, B->split);
        rel = B->Perr / (double) B->Pm[W - 1] + yerr / (double) y[W - 1];
        {
            slong top = (H[W] != 0) ? W : W - 1;
            int c;

            FLINT_ASSERT(H[top] != 0);
            c = flint_clz(H[top]);
            if (c != 0)
                mpn_lshift(H, H, top + 1, c);
            flint_mpn_copyi(B->Pm, H + (top + 1 - W), W);
            B->eP = B->eP + ey + FLINT_BITS * top - c;
            B->Perr = rev_new_err(B->Pm, W, rel, rev_ldexp(2.0 * W + 5.0, c));
        }
    }
#if FLINT_HAVE_FFT_SMALL
    B->QHok = 0;
#endif
}

/* ------------------------------------------------------------------------ */
/* the threaded region: the outputs of the current batch (jump) and the
   table pass of the next batch */

typedef struct
{
    rev_struct * B;
    int outputs;                /* the jump outputs of B->cur */
    int k;
    ulong n0;
    nn_ptr d[KMAX + 1];         /* (n0)_(2t) */
    slong dn[KMAX + 1];
    rev_pref_struct P0, Pk;        /* P_0 (input) and P_k (output) */
    rev_counter_struct next_out;
    const rev_batch_struct * T;    /* the table pass, or NULL */
    rev_counter_struct next_pow;
    rev_counter_struct next_sm;    /* the smooth powers (one task) */
}
rev_region_struct;

typedef struct
{
    rev_region_struct * R;
    slong wi;
}
rev_region_arg;

#if FLINT_HAVE_FFT_SMALL

/* P_t = P_0 y^t / d_t from the transform of P_0 */
static void
rev_jump_pref(rev_pref_struct * Pt, rev_region_struct * R, rev_worker_struct * w,
    int t, nn_ptr z, nn_ptr q)
{
    rev_struct * B = R->B;
    slong W = B->W, zn = 2 * W - B->flo, qn, drop;
    double rel;
    nn_srcptr yt = B->ym[t] + (B->Wy - W);

    rev_fft_fit_Y(B, t);
    fft_small_op_mul(w->X, B->QH, B->Y + t, B->FP);
    rev_fft_export(z, B, w->X);
    /* z is floor(P_0 y^t / B^flo) within 2 units; divide by d_t */
    flint_mpn_tdiv_q(q, z, zn, R->d[t], R->dn[t]);
    qn = rev_len(q, zn - R->dn[t] + 1);
    FLINT_ASSERT(qn >= W + 1);
    drop = qn - W;
    flint_mpn_copyi(Pt->m, q + drop, W);
    Pt->e = R->P0.e + B->ey[t] + FLINT_BITS * ((B->Wy - W) + B->flo + drop);
    rel = R->P0.err / (double) R->P0.m[W - 1] + (B->yerr[t] + 1.0) / (double) yt[W - 1];
    Pt->err = rev_new_err(Pt->m, W, rel, 4.0);
}

/* the output task c (the expensive numbers first: t = 1, ..., k - 1,
   then t = k, which is only the next prefactor, then t = 0) */
static void
rev_output_task(rev_region_struct * R, rev_worker_struct * w, slong c)
{
    rev_struct * B = R->B;
    const rev_batch_struct * T = &B->cur;
    slong W = B->W, wp = T->wd * B->db, s0 = wp % FLINT_BITS;
    slong off = W - 1 - B->flo;
    slong t;
    ulong n;
    nn_ptr scr, dl, z, q, Tm, U, pm;
    slong Tn, dn;
    double ed;
    rev_pref_struct P;

    t = (c < R->k - 1) ? c + 1 : (c == R->k - 1) ? R->k : 0;
    n = R->n0 - 2 * t;

    rev_wk_fit_ops(B, w);
    scr = rev_wk_tmp(w, rev_out_scratch(W, n, B->dmax));
    z = scr;                        /* 2 W limbs */
    q = z + 2 * W + 2;              /* 2 W + dmax limbs */
    pm = q + 2 * W + B->dmax + 4;   /* W */
    dl = pm + W + 2;                /* W + 2 */
    Tm = dl + W + 4;                /* W + 2 */
    U = Tm + W + 4;

    if (t == 0)
        P = R->P0;
    else
    {
        P.m = (t == R->k) ? R->Pk.m : pm;
        rev_jump_pref(&P, R, w, t, z, q);
        if (t == R->k)
        {
            /* the transform of the next P_0 */
            fft_small_fft_mpn(B->QHn, P.m, W, B->FP);
            R->Pk = P;
            return;
        }
    }

    rev_reduce(B, T, t);
    ed = rev_delta(dl, B->S[t], B->Sn[t], B->Serr[t], n, wp, W, U);
    dn = rev_len(dl, W);
    fft_small_fft_mpn(w->Z, dl, dn, B->FP);
    if (t == 0)
        fft_small_op_mul(w->Z, B->QH, w->Z, B->FP);
    else
    {
        fft_small_fft_mpn(w->X, P.m, W, B->FP);
        fft_small_op_mul(w->Z, w->X, w->Z, B->FP);
    }
    rev_fft_export_short(z, B, w->Z, dn);
    Tn = rev_T_from_window(Tm, z + off, W, s0);

    rev_round(B->num + t, B->den + t, n, &P, W, Tm, Tn, 2.0 * W + 4.0,
        ((double) P.m[W - 1] + 1.0) * ed, FLINT_BITS * (W - 1) - wp,
        B->den + t, U);
}

#endif

static void
rev_region_worker(void * varg)
{
    rev_region_arg * a = varg;
    rev_region_struct * R = a->R;
    rev_struct * B = R->B;
    rev_worker_struct * w = B->wk + a->wi;

#if FLINT_HAVE_FFT_SMALL
    if (R->outputs)
    {
        while (1)
        {
            slong c = rev_counter_fetch_add(&R->next_out, 1);
            if (c > R->k)
                break;
            rev_output_task(R, w, c);
        }
    }
#endif

    if (R->T != NULL)
    {
        const rev_batch_struct * T = R->T;
        int started = 0;

        if (B->sieve && rev_counter_fetch_add(&R->next_sm, 1) == 0)
            rev_smooth_pass(B, T);

        while (1)
        {
            slong i0 = rev_counter_fetch_add(&R->next_pow, TABLE_GRAB), i;
            if (i0 >= B->nr)
                break;
            if (!started)
            {
                rev_wk_table_start(B, w, T);
                started = 1;
            }
            if (B->mode == REV_MODE_GENERIC)
            {
                slong cnt = FLINT_MIN(i0 + TABLE_GRAB, B->nr) - i0;
                for (i = 0; i < cnt; i++)
                    rev_table_prep(B, B->ri[i0 + i], T, w->Aerr[T->par]);
                rev_kernel_generic_tiled(B, B->ri + i0, cnt, T->k, w->A[T->par],
                    T->sdlen, B->split, w->car);
                continue;
            }
            for (i = i0; i < FLINT_MIN(i0 + TABLE_GRAB, B->nr); i++)
            {
                rev_table_one(B, w, B->ri[i], T, NULL, 0, w->Aerr[T->par]);
                if (++w->cnt == CHUNK)
                    rev_wk_flush(B, w, T);
            }
        }
        if (started)
            rev_wk_flush(B, w, T);
    }
}

/* run the region on this thread and the workers obtained */
static void
rev_run_region(rev_t B, rev_region_struct * R, slong nworkers,
    thread_pool_handle * handles)
{
    slong i;
    int save;
    rev_region_arg * args;

    rev_fit_workers(B, nworkers + 1);
    args = flint_malloc((nworkers + 1) * sizeof(rev_region_arg));
    for (i = 0; i <= nworkers; i++)
    {
        args[i].R = R;
        args[i].wi = i;
    }
    /* no nested parallelism: the region uses the threads it has */
    save = flint_set_num_workers(0);
    for (i = 0; i < nworkers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, rev_region_worker, args + i + 1);
    rev_region_worker(args);
    flint_reset_num_workers(save);
    for (i = 0; i < nworkers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);
    flint_free(args);
}

static void
rev_region_init(rev_region_struct * R, rev_struct * B)
{
    memset(R, 0, sizeof(rev_region_struct));
    R->B = B;
    rev_counter_init(&R->next_out);
    rev_counter_init(&R->next_pow);
    rev_counter_init(&R->next_sm);
}

static void
rev_region_clear(rev_region_struct * R)
{
    rev_counter_clear(&R->next_out);
    rev_counter_clear(&R->next_pow);
    rev_counter_clear(&R->next_sm);
}

#if FLINT_HAVE_FFT_SMALL

/* the divisors (n0)_(2t), t <= k */
static nn_ptr
rev_divisors(rev_region_struct * R, slong dmax)
{
    nn_ptr dbuf = flint_malloc((R->k + 1) * (dmax + 2) * sizeof(ulong));
    ulong n0 = R->n0;
    int t;

    R->d[0] = NULL;
    R->dn[0] = 0;
    for (t = 1; t <= R->k; t++)
    {
        ulong a = n0 - 2 * (t - 1), hi, lo, ab[2];
        slong m = (t == 1) ? 1 : R->dn[t - 1];
        R->d[t] = dbuf + t * (dmax + 2);
        umul_ppmm(hi, lo, a, a - 1);
        ab[0] = lo;
        ab[1] = hi;
        if (t == 1)
        {
            R->d[t][0] = lo;
            R->d[t][1] = hi;
            R->dn[t] = 1 + (hi != 0);
        }
        else
        {
            slong abn = 1 + (hi != 0);
            if (m >= abn)
                flint_mpn_mul(R->d[t], R->d[t - 1], m, ab, abn);
            else
                flint_mpn_mul(R->d[t], ab, abn, R->d[t - 1], m);
            R->dn[t] = rev_len(R->d[t], m + abn);
        }
        FLINT_ASSERT(R->dn[t] <= dmax);
    }
    return dbuf;
}

#endif

/* ------------------------------------------------------------------------ */
/* the initial table */

/* The powers X_j = j^-nmax 2^wp (j = 2 i + 3, up to err[i] units) are
   built in limbs, entry by entry, without holding more than the table:
   the primes by powering in arb, the composites j = p q (p the smallest
   prime factor) as the products of the entries of p and q, truncated to
   the limbs that reach the frame.  The stages (the primes, then the
   composites in [L, 3 L) for L = 9, 27, ..., whose cofactors are below
   L) run on the threads given; with digits, a last stage converts the
   entries. */

typedef struct
{
    rev_struct * B;
    ulong nmax;
    slong wp;
    const unsigned int * spf;   /* the smallest prime factor of odd j */
    const slong * idx;          /* the entry of odd j, or -1 */
    const slong * todo;         /* the entries of the stage */
    slong ntodo;
    int stage;                  /* 0: primes, 1: composites, 2: digits */
}
rev_tab_struct;

/* a buffer for an entry of len limbs (in limbs or digits, with room for
   growth) */
static ulong *
rev_tab_alloc(rev_struct * B, slong len, slong * cap)
{
    *cap = (len * FLINT_BITS) / B->db + 1 + ROOM + 64;
    return flint_aligned_alloc(64, ((*cap * sizeof(ulong) + 63) / 64) * 64);
}

static void
rev_tab_prime(rev_tab_struct * A, slong i)
{
    rev_struct * B = A->B;
    ulong jj = B->jv[i];
    slong pj, len;
    nn_srcptr fd;
    arb_t u;
    fmpz_t f;
    int r1;

    arb_init(u);
    fmpz_init(f);
    pj = A->wp - A->nmax * rev_log2(jj) + 32;
    pj = FLINT_MAX(pj, 32);
    arb_ui_pow_ui(u, jj, A->nmax, pj);
    arb_inv(u, u, pj);
    r1 = arf_get_fmpz_fixed_si(f, arb_midref(u), -A->wp);
    len = fmpz_size(f);
    fd = COEFF_IS_MPZ(*f) ? COEFF_TO_PTR(*f)->_mp_d : (nn_srcptr) f;
    B->buf[i] = rev_tab_alloc(B, len, B->cap + i);
    flint_mpn_copyi(B->buf[i], fd, len);
    B->nd[i] = len;
    B->err[i] = rev_mag_2exp(arb_radref(u), A->wp) * 1.000001 + 2.0 + (r1 != 0);
    arb_clear(u);
    fmpz_clear(f);
}

/* X_j = floor(X_p X_q / 2^wp) from the top limbs of X_p and X_q */
static void
rev_tab_composite(rev_tab_struct * A, slong i)
{
    rev_struct * B = A->B;
    ulong jj = B->jv[i], p = A->spf[(jj - 3) / 2], q = jj / p;
    slong ip = A->idx[(p - 3) / 2], iq = A->idx[(q - 3) / 2], tl, len;
    nn_ptr z;

    FLINT_ASSERT(ip >= 0 && iq >= 0);
    tl = (slong) ((A->wp - A->nmax * rev_log2(jj)) / FLINT_BITS) + 1;
    z = flint_malloc((B->nd[ip] + B->nd[iq] + 1) * sizeof(ulong));
    B->err[i] = rev_mul_frame(z, &len, B->buf[ip], B->nd[ip], B->err[ip],
        B->buf[iq], B->nd[iq], B->err[iq], A->wp, tl);
    B->buf[i] = rev_tab_alloc(B, len, B->cap + i);
    flint_mpn_copyi(B->buf[i], z, len);
    B->nd[i] = len;
    flint_free(z);
}

#if REV_HAVE_IFMA
static void
rev_tab_digits(rev_tab_struct * A, slong i)
{
    rev_struct * B = A->B;
    slong len = B->nd[i];
    nn_ptr t;

    if (len == 0)
        return;
    t = flint_malloc(len * sizeof(ulong));
    flint_mpn_copyi(t, B->buf[i], len);
    B->nd[i] = rev_to_digits(B->buf[i], t, len, B->db);
    flint_free(t);
}
#endif

/* the entry todo[c] of the stage */
static void
rev_tab_task(slong c, void * arg)
{
    rev_tab_struct * A = arg;
    slong i = A->todo[c];

    if (A->stage == 0)
        rev_tab_prime(A, i);
    else if (A->stage == 1)
        rev_tab_composite(A, i);
#if REV_HAVE_IFMA
    else
        rev_tab_digits(A, i);
#endif
}

/* a stage on up to nthreads threads, the entries taken in order (the
   largest first) as the threads become free */
static void
rev_tab_run(rev_tab_struct * A, slong nthreads)
{
    flint_parallel_do(rev_tab_task, A, A->ntodo, FLINT_MAX(nthreads, 1),
        FLINT_PARALLEL_DYNAMIC);
}

static void
rev_init_table(rev_struct * B, ulong nmax, slong wp, ulong J, slong nthreads)
{
    slong no = (J >= 3) ? (J - 1) / 2 : 0, np, i, c;
    ulong L, p, m, j;
    unsigned int * spf;
    slong * todo, * idx;
    rev_tab_struct A;

    /* the entries: odd 3 <= j <= J, rough or smooth */
    B->J = J;
    B->jv = flint_malloc((no + 1) * sizeof(ulong));
    idx = flint_malloc((no + 1) * sizeof(slong));
    np = 0;
    for (i = 0; i < no; i++)
    {
        j = 2 * i + 3;
        if (rev_j_rough(j, B->sieve) || rev_j_smooth(j, B->sieve))
        {
            idx[i] = np;
            B->jv[np++] = j;
        }
        else
            idx[i] = -1;
    }
    B->np = np;
    B->buf = flint_malloc((np + 1) * sizeof(ulong *));
    B->off = flint_malloc(3 * (np + 1) * sizeof(slong));
    B->nd = B->off + np + 1;
    B->cap = B->nd + np + 1;
    B->err = flint_malloc((np + 1) * sizeof(double));
    B->ri = flint_malloc(2 * (np + 1) * sizeof(slong));
    B->si = B->ri + np + 1;
    B->nr = B->ns = 0;
    for (i = 0; i < np; i++)
    {
        if (rev_j_rough(B->jv[i], B->sieve))
            B->ri[B->nr++] = i;
        else
            B->si[B->ns++] = i;
    }

    if (np == 0)
    {
        flint_free(idx);
        return;
    }

    /* the smallest prime factors of the odd j <= J */
    spf = flint_calloc(no, sizeof(unsigned int));
    for (p = 3; p <= J; p += 2)
    {
        if (spf[(p - 3) / 2] != 0)
            continue;
        spf[(p - 3) / 2] = p;
        if (p <= J / p)
            for (m = p * p; m <= J; m += 2 * p)
                if (spf[(m - 3) / 2] == 0)
                    spf[(m - 3) / 2] = p;
    }

    todo = flint_malloc(np * sizeof(slong));
    A.B = B;
    A.nmax = nmax;
    A.wp = wp;
    A.spf = spf;
    A.idx = idx;
    A.todo = todo;

    /* the primes, largest entries (smallest j) first */
    for (i = c = 0; i < np; i++)
        if (spf[(B->jv[i] - 3) / 2] == B->jv[i])
            todo[c++] = i;
    A.ntodo = c;
    A.stage = 0;
    rev_tab_run(&A, nthreads);

    /* the composites, by ranges [L, 3 L) (the factors of an entry being
       entries) */
    for (L = 9; L <= J; L *= 3)
    {
        c = 0;
        for (i = 0; i < np; i++)
            if (B->jv[i] >= L && B->jv[i] < 3 * L && spf[(B->jv[i] - 3) / 2] != B->jv[i])
                todo[c++] = i;
        A.ntodo = c;
        A.stage = 1;
        rev_tab_run(&A, nthreads);
        if (L > J / 3)
            break;
    }

#if REV_HAVE_IFMA
    /* the rough powers to digits (the smooth ones stay in limbs) */
    if (B->mode == REV_MODE_IFMA)
    {
        for (i = 0; i < B->nr; i++)
            todo[i] = B->ri[i];
        A.ntodo = B->nr;
        A.stage = 2;
        rev_tab_run(&A, nthreads);
    }
#endif

    for (i = 0; i < np; i++)
        B->off[i] = 0;

    flint_free(todo);
    flint_free(spf);
    flint_free(idx);
}

/* ------------------------------------------------------------------------ */
/* init and clear */

static void
rev_init(rev_t B, ulong nmax, int flags)
{
    slong wp, W, prec, t, i;
    ulong J;
    arb_t u, v;

    memset(B, 0, sizeof(rev_struct));
    for (t = 0; t < KMAX; t++)
    {
        fmpz_init(B->num + t);
        fmpz_init(B->den + t);
    }

    B->split = (flags & BERNOULLI_REV_TEST_SPLIT) != 0;

    /* the kernel, from the largest j at the precision of the generic
       frame (the digit frames are at most a limb larger) */
    wp = bernoulli_global_prec(nmax) + GUARD + FLINT_BITS;
    J = bernoulli_zeta_terms(nmax - 2 * (rev_batch_k(nmax) - 1), wp);
    B->mode = REV_MODE_GENERIC;
#if REV_HAVE_IFMA
    if (J < (UWORD(1) << 24) && !(flags & BERNOULLI_REV_TEST_GENERIC))
        B->mode = REV_MODE_IFMA;
#endif
    B->db = (B->mode == REV_MODE_IFMA) ? 50 : FLINT_BITS;

    /* the sieve: forced by the flags, else for large n */
    if (flags & BERNOULLI_REV_SIEVE_15)
        B->sieve = 2;
    else if (flags & BERNOULLI_REV_SIEVE_3)
        B->sieve = 1;
    else if (flags & BERNOULLI_REV_NO_SIEVE)
        B->sieve = 0;
    else
        B->sieve = (nmax >= SIEVE_CUTOFF) ? 2 : 0;

    B->wd = rev_frame_digits(B, nmax);
    wp = B->wd * B->db;
    J = bernoulli_zeta_terms(nmax - 2 * (rev_batch_k(nmax) - 1), wp);
    arb_init(u);
    arb_init(v);

    rev_init_table(B, nmax, wp, J, (rev_prod_limbs(wp) >= THREAD_CUTOFF)
        ? flint_get_num_threads() : 1);

    B->Salloc = wp / FLINT_BITS + ROOM + 64;
    B->Sdalloc = B->wd + ROOM + 64;
    for (t = 0; t < KMAX; t++)
        B->S[t] = flint_malloc(B->Salloc * sizeof(ulong));
    if (B->sieve)
        for (t = 0; t < KMAX; t++)
            for (i = 0; i < 2; i++)
                B->G[i][t] = flint_malloc(B->Salloc * sizeof(ulong));
    rev_fit_workers(B, 1);

    /* divisors (n0)_(2t), t <= KMAX */
    B->dmax = (slong) (2 * KMAX * rev_log2(nmax + 1.0) / FLINT_BITS) + 3;

    /* the prefactor 2 nmax! / (2 pi)^nmax and the powers (2 pi)^(2t) */
    W = rev_prod_limbs(wp);
    B->W = W;
    B->Wy = W;
    prec = FLINT_BITS * W + 64;
    B->Pm = flint_malloc((W + 2) * sizeof(ulong));
    arb_const_pi(u, prec);
    arb_mul_2exp_si(u, u, 1);
    arb_mul(u, u, u, prec);
    arb_set(v, u);
    for (t = 1; t <= KMAX; t++)
    {
        B->ym[t] = flint_malloc((W + 2) * sizeof(ulong));
        rev_arb_get_mpn(B->ym[t], B->ey + t, B->yerr + t, v, W);
        arb_mul(v, v, u, prec);
    }
    arb_pow_ui(v, u, nmax / 2, prec);
    arb_fac_ui(u, nmax, prec);
    arb_mul_2exp_si(u, u, 1);
    arb_div(u, u, v, prec);
    rev_arb_get_mpn(B->Pm, &B->eP, &B->Perr, u, W);

#if FLINT_HAVE_FFT_SMALL
    rev_fft_setup(B);
#endif

    arb_clear(u);
    arb_clear(v);
}

static void
rev_clear(rev_t B)
{
    slong i, t;
    for (t = 0; t < KMAX; t++)
    {
        fmpz_clear(B->num + t);
        fmpz_clear(B->den + t);
    }
    for (i = 0; i < B->np; i++)
        flint_aligned_free(B->buf[i]);
    flint_free(B->buf);
    flint_free(B->off);
    flint_free(B->err);
    flint_free(B->jv);
    flint_free(B->ri);
    if (B->sieve)
        for (t = 0; t < KMAX; t++)
        {
            flint_free(B->G[0][t]);
            flint_free(B->G[1][t]);
        }
    for (t = 0; t < KMAX; t++)
        flint_free(B->S[t]);
    for (t = 1; t <= KMAX; t++)
        flint_free(B->ym[t]);
#if FLINT_HAVE_FFT_SMALL
    rev_fft_clear(B);
#endif
    rev_clear_workers(B);
    flint_free(B->Pm);
}

/* ------------------------------------------------------------------------ */
/* one batch */

static void
rev_refill(rev_t B, ulong n)
{
    slong W_new, d, t, nworkers = 0, wp;
    thread_pool_handle * handles = NULL;
    rev_region_struct R;
    ulong n1;
    int threaded;

    /* the table pass of this batch, unless it was done with the previous
       outputs */
    if (B->have_next)
        B->cur = B->next;
    else
    {
        rev_prepare(B, &B->cur, n);
        rev_region_init(&R, B);
        R.T = &B->cur;
        for (t = 0; t < B->nwk; t++)
            B->wk[t].touched[B->cur.par] = 0;
        if (B->W >= THREAD_CUTOFF)
            nworkers = flint_request_threads(&handles, flint_get_num_threads());
        rev_run_region(B, &R, nworkers, handles);
        flint_give_back_threads(handles, nworkers);
        rev_region_clear(&R);
        nworkers = 0;
        handles = NULL;
    }
    wp = B->cur.wd * B->db;

    /* lower the product precision in steps of about 3% */
    W_new = rev_prod_limbs(wp);
    if (W_new < B->W - FLINT_MAX(2, B->W / 32))
    {
        d = B->W - W_new;
        flint_mpn_copyi(B->Pm, B->Pm + d, W_new);
        B->eP += FLINT_BITS * d;
        B->Perr = rev_ldexp(B->Perr, -FLINT_BITS * d) * 1.000001 + 1.0;
        B->W = W_new;
#if FLINT_HAVE_FFT_SMALL
        rev_fft_setup(B);
#endif
    }

    /* the denominators */
    for (t = 0; t < B->cur.k; t++)
        arith_bernoulli_number_denom(B->den + t, B->cur.n0 - 2 * t);

    /* prepare the next table pass */
    n1 = B->cur.n0 - 2 * B->cur.k;
    B->have_next = (n1 >= BERNOULLI_REV_MIN);
    if (B->have_next)
    {
        rev_prepare(B, &B->next, n1);
        for (t = 0; t < B->nwk; t++)
            B->wk[t].touched[B->next.par] = 0;
    }

    /* threads, if any are available now */
#if FLINT_HAVE_FFT_SMALL
    if (B->fft && B->W >= THREAD_CUTOFF)
        nworkers = flint_request_threads(&handles, flint_get_num_threads());
#endif
    threaded = (nworkers > 0);

    rev_region_init(&R, B);
    R.T = B->have_next ? &B->next : NULL;

    if (!threaded)
    {
        int save = flint_set_num_workers(0);
        rev_outputs_chain(B, B->wk);
        flint_reset_num_workers(save);
        /* the table pass on the threads available now (without
           fft_small, or when none were free for the outputs) */
        if (R.T != NULL)
        {
            flint_give_back_threads(handles, nworkers);
            handles = NULL;
            nworkers = 0;
            if (B->W >= THREAD_CUTOFF)
                nworkers = flint_request_threads(&handles, flint_get_num_threads());
            rev_run_region(B, &R, nworkers, handles);
        }
    }
#if FLINT_HAVE_FFT_SMALL
    else
    {
        nn_ptr dbuf, pk;

        R.outputs = 1;
        R.k = B->cur.k;
        R.n0 = B->cur.n0;
        dbuf = rev_divisors(&R, B->dmax);
        R.P0.m = B->Pm;
        R.P0.e = B->eP;
        R.P0.err = B->Perr;
        pk = flint_malloc((B->W + 2) * sizeof(ulong));
        R.Pk.m = pk;
        if (!B->QHok)
            fft_small_fft_mpn(B->QH, B->Pm, B->W, B->FP);

        rev_run_region(B, &R, nworkers, handles);

        FLINT_SWAP(fft_small_op_struct, *B->QH, *B->QHn);
        B->QHok = 1;
        flint_mpn_copyi(B->Pm, R.Pk.m, B->W);
        B->eP = R.Pk.e;
        B->Perr = R.Pk.err;
        flint_free(pk);
        flint_free(dbuf);
    }
#endif
    flint_give_back_threads(handles, nworkers);
    rev_region_clear(&R);

    B->nbuf = B->cur.k;
    B->bpos = 0;
}

/* ------------------------------------------------------------------------ */
/* the interface */

void
_bernoulli_rev_init_flags(bernoulli_rev_t iter, ulong nmax, int flags)
{
    if (nmax > BERNOULLI_REV_NMAX)
        flint_throw(FLINT_ERROR, "bernoulli_rev_init: n = %wu is too large\n", nmax);

    nmax -= nmax % 2;
    iter->n = nmax;
    iter->data = NULL;

    if (nmax >= BERNOULLI_REV_MIN)
    {
        iter->data = flint_malloc(sizeof(rev_struct));
        rev_init(iter->data, nmax, flags);
    }
}

void
bernoulli_rev_init(bernoulli_rev_t iter, ulong nmax)
{
    _bernoulli_rev_init_flags(iter, nmax, 0);
}

void
bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernoulli_rev_t iter)
{
    ulong n = iter->n;
    rev_struct * B = iter->data;

    if (n < BERNOULLI_REV_MIN)
    {
        _bernoulli_fmpq_ui(numer, denom, n);
        if (n != 0)
            iter->n -= 2;
        return;
    }

    if (B->bpos == B->nbuf)
        rev_refill(B, n);

    fmpz_swap(numer, B->num + B->bpos);
    fmpz_swap(denom, B->den + B->bpos);
    B->bpos++;
    iter->n -= 2;
}

void
bernoulli_rev_clear(bernoulli_rev_t iter)
{
    if (iter->data != NULL)
    {
        rev_clear(iter->data);
        flint_free(iter->data);
        iter->data = NULL;
    }
}
