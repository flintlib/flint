/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "thread_support.h"
#include "impl.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/*
    Transform-space versions of the matrix operations of the half-gcd.

    The 2 x 2 products have four outputs, each a sum of two products, so
    with fft_small every operand is transformed once and every output is
    transformed back once:

      matrix_mul (M <- M M1)        8 forward, 4 inverse transforms
                                    (Strassen-like: 7 products, 21)
      matrix_adjust (low parts)     6 forward, 2 inverse (4 products: 12)
      matrix_adjust_mul (both of    10 forward, 6 inverse, the entries of
        the above, sharing M1)      M1 transformed once (separately: 20)
      mul_matrix_vector (gcdext)    6 forward, 2 inverse (4 products: 12)
      matrix_apply (mod B^N - 1)    7 forward, 2 inverse cyclic transforms
                                    of half length (4 cyclic products: 12)

    The entries are nonnegative. Outputs that are sums of two products are
    exported unsigned; the differences of matrix_adjust are exported
    signed; the differences of matrix_apply are made nonnegative modulo
    B^N - 1 by replacing -x with (B^N - 1) - x, whose transform is the
    transform of B^N - 1 (all ones) minus that of x.

    Each operation runs in two stages, each a set of independent jobs: the
    forward transforms, and for each output (or pair of outputs) its
    pointwise combination, inverse transform and export. For large
    operands, the jobs of a stage run on several threads; fft_small also
    threads within each transform from some size on, but not in the
    two-prime geometry used up to FLINT_MPN_MUL_NP2_MAX_AN limbs, and over
    at most its number of primes. The forward transforms are listed with
    the longer and shorter operands interleaved, so that the strided
    distribution of the jobs over the threads is balanced.
*/

/* Whether products of lengths at least an >= bn should go through
   fft_small when the transforms are shared as above: the rule of
   flint_mpn_mul with the threshold T replaced by the given cutoff. */
static int
use_fft(mp_size_t an, mp_size_t bn, mp_size_t cutoff)
{
    if (an < bn)
        FLINT_SWAP(mp_size_t, an, bn);
    return (bn >= cutoff) || (an + bn > 2 * cutoff && 3 * bn > cutoff);
}

/* The cutoffs do not follow FLINT_FFT_SMALL_MUL_THRESHOLD: tune-gcd gives
   the same optimum, a flat one around 150 limbs, on x86 machines whose
   thresholds are 240 and 500. Platforms with a higher threshold (generic,
   arm64), which have not been measured, keep 0.3 times it. */
#define FLINT_MPN_HGCD_FFT_CUTOFF_DEFAULT \
    FLINT_MAX(150, 3 * FLINT_FFT_SMALL_MUL_THRESHOLD / 10)
#ifndef FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF
# define FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF FLINT_MPN_HGCD_FFT_CUTOFF_DEFAULT
#endif
#ifndef FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF
# define FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF FLINT_MPN_HGCD_FFT_CUTOFF_DEFAULT
#endif
#ifndef FLINT_MPN_HGCD_APPLY_FFT_CUTOFF
# define FLINT_MPN_HGCD_APPLY_FFT_CUTOFF FLINT_MPN_HGCD_FFT_CUTOFF_DEFAULT
#endif

/* Whether count operands of the plan fit within the global bound on
   transform storage; if not, the functions below decline, and their
   callers fall back to separate products, which need less memory. */
static int
pool_fits(const fft_small_plan_struct * P, ulong count)
{
    return count * fft_small_op_sizeof_data(P) <= flint_fft_small_max_transformed_ring_size;
}

/* z (zn limbs) = the exact nonnegative value of Z, a sum of at most terms
   products of operands whose products have at most nslots slots; the
   value must fit in zn limbs. */
static void
export_unsigned(mp_ptr z, mp_size_t zn, fft_small_op_struct * Z,
    ulong terms, ulong nslots, const fft_small_plan_struct * P)
{
    ulong need = fft_small_export_mpn_limbs(P, terms, nslots);
    ulong en = FLINT_MIN((ulong) zn, need);

    fft_small_ifft_trunc(Z, nslots, n_cdiv(en * FLINT_BITS, P->bits), NULL, P);
    fft_small_export_mpn(z, en, Z, P);
    if ((ulong) zn > en)
        flint_mpn_zero(z + en, zn - en);
}

/* (buf, returned length) = |Z|, sign in *sgn; buf needs
   fft_small_export_mpn_signed_limbs(P, nslots) limbs. */
static mp_size_t
export_signed(mp_ptr buf, int * sgn, fft_small_op_struct * Z,
    ulong nslots, const fft_small_plan_struct * P)
{
    ulong need = fft_small_export_mpn_signed_limbs(P, nslots);
    mp_size_t len;

    fft_small_ifft_trunc(Z, nslots, 0, NULL, P);
    fft_small_export_mpn_signed(buf, need, sgn, Z, nslots, P);

    len = need;
    MPN_NORM(buf, len);
    return len;
}

/* Whether to run the jobs of each stage on several threads, for
   operands of len limbs. The cutoff depends mostly on the cost of waking
   threads: 250 is optimal on a Zen 3 machine, 1000 on a virtual machine
   where an empty flint_parallel_do takes 27 us. */
#ifndef FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF
# define FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF 500
#endif

static int
use_parallel(mp_size_t len)
{
    return len >= FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF && flint_get_num_threads() > 1;
}

static void
run_jobs(do_func_t f, void * args, slong n, int par)
{
    slong i;

    if (par)
        flint_parallel_do(f, args, n, 0, FLINT_PARALLEL_STRIDED);
    else
        for (i = 0; i < n; i++)
            f(i, args);
}

/* Forward transforms: of (a, an); or of a folded to N limbs mod B^N - 1
   in tmp (fold != 0); or of B^N - 1 written to tmp (ones != 0). */
typedef struct
{
    fft_small_op_struct * X;
    mp_srcptr a;
    mp_size_t an;
    mp_ptr tmp;
    mp_size_t N;
    int fold;
    int ones;
}
fwd_job_struct;

typedef struct
{
    fwd_job_struct * jobs;
    const fft_small_plan_struct * P;
}
fwd_args_struct;

static void fold_bnm1(mp_ptr r, mp_srcptr v, mp_size_t vn, mp_size_t N);

static void
fwd_worker(slong i, void * varg)
{
    fwd_args_struct * A = (fwd_args_struct *) varg;
    fwd_job_struct * J = A->jobs + i;

    if (J->fold)
    {
        fold_bnm1(J->tmp, J->a, J->an, J->N);
        fft_small_fft_mpn(J->X, J->tmp, J->N, A->P);
    }
    else if (J->ones)
    {
        mp_size_t k;
        for (k = 0; k < J->N; k++)
            J->tmp[k] = UWORD_MAX;
        fft_small_fft_mpn(J->X, J->tmp, J->N, A->P);
    }
    else
    {
        fft_small_fft_mpn(J->X, J->a, J->an, A->P);
    }
}

static void
fwd_set(fwd_job_struct * J, fft_small_op_struct * X, mp_srcptr a, mp_size_t an)
{
    J->X = X;
    J->a = a;
    J->an = an;
    J->tmp = NULL;
    J->N = 0;
    J->fold = 0;
    J->ones = 0;
}

static void
forward(fwd_job_struct * jobs, slong n, const fft_small_plan_struct * P, int par)
{
    fwd_args_struct A;
    A.jobs = jobs;
    A.P = P;
    run_jobs(fwd_worker, &A, n, par);
}

/* Pointwise stage: each job evaluates one or two sums
   Z = X0 Y0 +- X1 Y1 +- X2 Y2 in order; Z may be X0 but no other operand
   of its sum. */
typedef struct
{
    fft_small_op_struct * Z;
    const fft_small_op_struct * X[3];
    const fft_small_op_struct * Y[3];
    int neg[3];
    int nterms;
}
dot_struct;

typedef struct
{
    dot_struct d[2];
    int nd;
}
dot_job_struct;


static void
dot_set2(dot_struct * d, fft_small_op_struct * Z,
    const fft_small_op_struct * X0, const fft_small_op_struct * Y0,
    int neg1, const fft_small_op_struct * X1, const fft_small_op_struct * Y1)
{
    d->Z = Z;
    d->X[0] = X0; d->Y[0] = Y0; d->neg[0] = 0;
    d->X[1] = X1; d->Y[1] = Y1; d->neg[1] = neg1;
    d->nterms = 2;
}

static void
dot_run(const dot_job_struct * J, const fft_small_plan_struct * P)
{
    int k, t;

    for (k = 0; k < J->nd; k++)
    {
        const dot_struct * d = J->d + k;

        fft_small_op_mul(d->Z, d->X[0], d->Y[0], P);
        for (t = 1; t < d->nterms; t++)
        {
            if (d->neg[t])
                fft_small_op_submul(d->Z, d->X[t], d->Y[t], P);
            else
                fft_small_op_addmul(d->Z, d->X[t], d->Y[t], P);
        }
    }
}


/* Export stage. */
#define EXPORT_UNSIGNED 0   /* z (zn limbs) = Z, a sum of terms products */
#define EXPORT_ADJUST 1     /* x <- x_hi B^p + Z (signed), see adjust_part */
#define EXPORT_CYCLIC 2     /* x (nn limbs) <- Z mod B^N - 1, see below */

typedef struct
{
    int kind;
    fft_small_op_struct * Z;
    mp_ptr z;
    mp_size_t zn;
    ulong nslots;
    ulong terms;
    /* EXPORT_ADJUST */
    mp_size_t n;
    mp_size_t p;
    mp_size_t mn;
    mp_ptr buf;
    mp_limb_t hi;       /* output: the carry out of the top limb */
    /* EXPORT_CYCLIC */
    mp_size_t N;
    mp_size_t nn;
    mp_ptr t;           /* N limbs, and buf of zn limbs */
}
exp_job_struct;

/* Output stage: each job computes the pointwise products of its outputs
   (dots), then exports them (ex). */
typedef struct
{
    dot_job_struct dots;
    exp_job_struct ex[2];
    int nex;
}
out_job_struct;

typedef struct
{
    out_job_struct * jobs;
    const fft_small_plan_struct * P;
}
out_args_struct;

/* x = x_hi B^p + S, the final stage of the matrix adjustment for one of
   a and b: S (signed, of at most p + mn limbs) is exported from Z, whose
   products have at most nslots slots, into buf (which needs
   fft_small_export_mpn_signed_limbs(P, nslots) limbs). x has n limbs.
   Returns the carry out of limb n - 1. */
static mp_limb_t
adjust_part(fft_small_op_struct * Z, ulong nslots, mp_size_t n, mp_ptr xp,
    mp_size_t p, mp_size_t mn, mp_ptr buf, const fft_small_plan_struct * P)
{
    mp_size_t sn;
    int sgn;

    sn = export_signed(buf, &sgn, Z, nslots, P);
    FLINT_ASSERT(sn <= p + mn);
    (void) mn;
    flint_mpn_zero(xp, p);

    if (sn == 0)
        return 0;
    else if (!sgn)
        return mpn_add(xp, xp, n, buf, sn);
    else
    {
        mp_limb_t cy = mpn_sub(xp, xp, n, buf, sn);
        FLINT_ASSERT(cy == 0);
        (void) cy;
        return 0;
    }
}

/* the new size of (a, b) after adjust_part with carries ah, bh */
static mp_size_t
adjust_size(mp_size_t n, mp_ptr ap, mp_ptr bp, mp_limb_t ah, mp_limb_t bh)
{
    if (ah > 0 || bh > 0)
    {
        ap[n] = ah;
        bp[n] = bh;
        n++;
    }
    else
    {
        /* The subtraction can reduce the size by at most one limb. */
        if (ap[n - 1] == 0 && bp[n - 1] == 0)
            n--;
    }
    FLINT_ASSERT(ap[n - 1] > 0 || bp[n - 1] > 0);
    return n;
}

static void
exp_run(exp_job_struct * J, const fft_small_plan_struct * P)
{
    if (J->kind == EXPORT_UNSIGNED)
    {
        export_unsigned(J->z, J->zn, J->Z, J->terms, J->nslots, P);
    }
    else if (J->kind == EXPORT_ADJUST)
    {
        J->hi = adjust_part(J->Z, J->nslots, J->n, J->z, J->p, J->mn, J->buf, P);
    }
    else
    {
        /* the cyclic product, reduced mod B^N - 1, is known to fit in nn
           limbs */
        fft_small_ifft(J->Z, P);
        fft_small_export_mpn(J->buf, J->zn, J->Z, P);
        fold_bnm1(J->t, J->buf, J->zn, J->N);
        FLINT_ASSERT(flint_mpn_zero_p(J->t + J->nn, J->N - J->nn));
        flint_mpn_copyi(J->z, J->t, J->nn);
    }
}

static void
out_worker(slong i, void * varg)
{
    out_args_struct * A = (out_args_struct *) varg;
    out_job_struct * J = A->jobs + i;
    int k;

    dot_run(&J->dots, A->P);
    for (k = 0; k < J->nex; k++)
        exp_run(J->ex + k, A->P);
}

static void
outputs(out_job_struct * jobs, slong n, const fft_small_plan_struct * P, int par)
{
    out_args_struct A;
    A.jobs = jobs;
    A.P = P;
    run_jobs(out_worker, &A, n, par);
}

static void
exp_set_unsigned(exp_job_struct * J, mp_ptr z, mp_size_t zn,
    fft_small_op_struct * Z, ulong terms, ulong nslots)
{
    J->kind = EXPORT_UNSIGNED;
    J->z = z;
    J->zn = zn;
    J->Z = Z;
    J->terms = terms;
    J->nslots = nslots;
}

static void
exp_set_adjust(exp_job_struct * J, mp_ptr xp, fft_small_op_struct * Z,
    ulong nslots, mp_size_t n, mp_size_t p, mp_size_t mn, mp_ptr buf)
{
    J->kind = EXPORT_ADJUST;
    J->z = xp;
    J->Z = Z;
    J->nslots = nslots;
    J->n = n;
    J->p = p;
    J->mn = mn;
    J->buf = buf;
    J->hi = 0;
}

/* The size of M after the four new entries of at most rn + mn + 1 limbs
   have been written. */
static void
matrix_set_size(flint_mpn_hgcd_matrix_struct * M, mp_size_t rn, mp_size_t mn)
{
    mp_size_t n;

    /* Index of last potentially non-zero limb, size is one greater. For
       consecutive half-gcd matrices this loop runs at most three times
       (see _flint_mpn_hgcd_matrix_mul), but other products may shrink
       more. */
    n = rn + mn;
    while (n > 0 && (M->p[0][0][n] | M->p[0][1][n] | M->p[1][0][n] | M->p[1][1][n]) == 0)
        n--;

    FLINT_ASSERT((M->p[0][0][n] | M->p[0][1][n] | M->p[1][0][n] | M->p[1][1][n]) > 0);

    M->n = n + 1;
}

/* The forward transforms of the entries of M (rn limbs) into X and of M1
   (mn limbs) into Y, interleaved: 8 jobs. */
static void
fwd_set_matrices(fwd_job_struct * jobs,
    fft_small_op_struct X[2][2], const flint_mpn_hgcd_matrix_struct * M, mp_size_t rn,
    fft_small_op_struct Y[2][2], const flint_mpn_hgcd_matrix_struct * M1, mp_size_t mn)
{
    int i, j, k = 0;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fwd_set(jobs + k++, X[i] + j, M->p[i][j], rn);
            fwd_set(jobs + k++, Y[i] + j, M1->p[i][j], mn);
        }
}

/* The output stage of M <- M M1, one job per row of M, with the
   temporaries T0, T1:

        r00' = r00 m00 + r01 m10  -> T0
        r01' = r00 m01 + r01 m11  -> X00 (r00, r01 dead)
        r10' = r10 m00 + r11 m10  -> T1
        r11' = r10 m01 + r11 m11  -> X10 (r10, r11 dead)

   The entries of the product have at most rn + mn + 1 limbs; the
   allocation of M allows this. */
static void
out_set_matmul(out_job_struct * jobs, flint_mpn_hgcd_matrix_struct * M,
    fft_small_op_struct X[2][2], fft_small_op_struct Y[2][2],
    fft_small_op_struct * T0, fft_small_op_struct * T1,
    ulong nslots, mp_size_t rn, mp_size_t mn)
{
    mp_size_t zn = rn + mn + 1;
    fft_small_op_struct * T[2];
    int r;

    FLINT_ASSERT(zn <= M->alloc);
    T[0] = T0;
    T[1] = T1;

    for (r = 0; r < 2; r++)
    {
        out_job_struct * J = jobs + r;

        dot_set2(J->dots.d + 0, T[r], X[r] + 0, Y[0] + 0, 0, X[r] + 1, Y[1] + 0);
        dot_set2(J->dots.d + 1, X[r] + 0, X[r] + 0, Y[0] + 1, 0, X[r] + 1, Y[1] + 1);
        J->dots.nd = 2;
        exp_set_unsigned(J->ex + 0, M->p[r][0], zn, T[r], 2, nslots);
        exp_set_unsigned(J->ex + 1, M->p[r][1], zn, X[r] + 0, 2, nslots);
        J->nex = 2;
    }
}

/* one output of the adjustment: Z = X0 Y0 - X1 Y1, added to x */
static void
out_set_adjust(out_job_struct * J, mp_ptr xp, fft_small_op_struct * Z,
    const fft_small_op_struct * X0, const fft_small_op_struct * Y0,
    const fft_small_op_struct * X1, const fft_small_op_struct * Y1,
    ulong nslots, mp_size_t n, mp_size_t p, mp_size_t mn, mp_ptr buf)
{
    dot_set2(J->dots.d, Z, X0, Y0, 1, X1, Y1);
    J->dots.nd = 1;
    exp_set_adjust(J->ex, xp, Z, nslots, n, p, mn, buf);
    J->nex = 1;
}

/******************************************************************************

    M <- M * M1

******************************************************************************/

int
_flint_mpn_hgcd_matrix_mul_fft(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], Y[2][2], T[2];
    fwd_job_struct fj[8];
    out_job_struct oj[2];
    mp_size_t rn = M->n, mn = M1->n;
    ulong nslots;
    int i, j, par;

    if (!use_fft(rn, mn, FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF))
        return 0;

    R = get_default_mpn_ctx();

    /* Each output is a sum of two products of nonnegative entries. */
    if (!fft_small_plan_init_mpn(P, R, rn, mn, 2, 0))
        return 0;

    if (!pool_fits(P, 10))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots = fft_small_plan_mpn_nslots(P, rn, mn);
    par = use_parallel(FLINT_MAX(rn, mn));

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 10, 0,
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_op_init_borrowed(Y[i] + j, P, fft_small_op_pool_take(pool));
        }
        fft_small_op_init_borrowed(T + i, P, fft_small_op_pool_take(pool));
    }

    fwd_set_matrices(fj, X, M, rn, Y, M1, mn);
    forward(fj, 8, P, par);

    out_set_matmul(oj, M, X, Y, T + 0, T + 1, nslots, rn, mn);
    outputs(oj, 2, P, par);
    matrix_set_size(M, rn, mn);

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            fft_small_op_pool_give(pool, X[i][j].data);
            fft_small_op_pool_give(pool, Y[i][j].data);
        }
        fft_small_op_pool_give(pool, T[i].data);
    }
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);
    return 1;
}

/******************************************************************************

    (a; b) <- M^{-1} (a; b) on the low p limbs

******************************************************************************/

/*
    Sets a = a_hi B^p + (m11 a_lo - m01 b_lo), b = b_hi B^p + (m00 b_lo -
    m10 a_lo) where a_lo, b_lo are the low p limbs; see
    _flint_mpn_hgcd_matrix_adjust. Returns the new size, or 0 if fft_small
    is not used.
*/
mp_size_t
_flint_mpn_hgcd_matrix_adjust_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], A, B;
    fwd_job_struct fj[6];
    out_job_struct oj[2];
    mp_size_t mn = M->n, an, bn;
    mp_ptr buf;
    ulong nslots, need;
    int i, j, par;

    if (!use_fft(p, mn, FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF))
        return 0;

    an = bn = p;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);
    if (an == 0 || bn == 0)
        return 0;

    R = get_default_mpn_ctx();

    /* Each output is a difference of two products. */
    if (!fft_small_plan_init_mpn(P, R, p, mn, 2, 1))
        return 0;

    if (!pool_fits(P, 6))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots = fft_small_plan_mpn_nslots(P, p, mn);
    need = fft_small_export_mpn_signed_limbs(P, nslots);
    par = use_parallel(FLINT_MAX(p, mn));

    /* the export buffers (one for each of a and b) are the pool's staging
       region, which like the slabs comes from the thread's fft_small
       buffer */
    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 6,
                          2 * need * sizeof(mp_limb_t),
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));

    fwd_set(fj + 0, &A, ap, an);
    fwd_set(fj + 1, X[0] + 0, M->p[0][0], mn);
    fwd_set(fj + 2, &B, bp, bn);
    fwd_set(fj + 3, X[0] + 1, M->p[0][1], mn);
    fwd_set(fj + 4, X[1] + 0, M->p[1][0], mn);
    fwd_set(fj + 5, X[1] + 1, M->p[1][1], mn);
    forward(fj, 6, P, par);

    /* a' = m11 a - m01 b -> X11, b' = m00 b - m10 a -> X00 */
    buf = fft_small_op_pool_scratch_take(pool, 2 * need, 0);
    out_set_adjust(oj + 0, ap, X[1] + 1, X[1] + 1, &A, X[0] + 1, &B, nslots, n, p, mn, buf);
    out_set_adjust(oj + 1, bp, X[0] + 0, X[0] + 0, &B, X[1] + 0, &A, nslots, n, p, mn, buf + need);
    outputs(oj, 2, P, par);
    n = adjust_size(n, ap, bp, oj[0].ex[0].hi, oj[1].ex[0].hi);
    fft_small_op_pool_scratch_give(pool, buf);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_pool_give(pool, X[i][j].data);
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    return n;
}

/******************************************************************************

    (a; b) <- M1^{-1} (a; b) on the low p limbs, then M <- M * M1

******************************************************************************/

/*
    The two operations following the second recursive call of the half-gcd,
    which share the matrix M1: with entries of M of about the size of a_lo
    and b_lo, one plan serves both, and the entries of M1 are transformed
    once (10 forward transforms instead of 14, out of 20 in total).
    Returns the new size of a and b, or 0 (doing nothing) if fft_small is
    not used, in which case the caller performs the two operations
    separately.
*/
mp_size_t
_flint_mpn_hgcd_matrix_adjust_mul_fft(flint_mpn_hgcd_matrix_struct * M,
    const flint_mpn_hgcd_matrix_struct * M1,
    mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], Y[2][2], A, B, Za, Zb, T[2];
    fwd_job_struct fj[10];
    out_job_struct oj[4];
    mp_size_t rn = M->n, mn = M1->n, an, bn, lmax, lmin;
    mp_ptr buf;
    ulong nslots_adj, nslots_mul, need;
    int i, j, par;

    if (!use_fft(p, mn, FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF) ||
        !use_fft(rn, mn, FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF))
        return 0;

    /* one transform length for both products: only when it does not
       penalize the shorter one much */
    lmax = FLINT_MAX(p, rn) + mn;
    lmin = FLINT_MIN(p, rn) + mn;
    if (4 * lmax > 5 * lmin)
        return 0;

    an = bn = p;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);
    if (an == 0 || bn == 0)
        return 0;

    R = get_default_mpn_ctx();

    /* differences of two products for the adjustment, sums of two products
       for the matrix product */
    if (!fft_small_plan_init_mpn(P, R, FLINT_MAX(p, rn), mn, 2, 1))
        return 0;

    if (!pool_fits(P, 14))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots_adj = fft_small_plan_mpn_nslots(P, p, mn);
    nslots_mul = fft_small_plan_mpn_nslots(P, rn, mn);
    need = fft_small_export_mpn_signed_limbs(P, nslots_adj);
    par = use_parallel(FLINT_MAX(FLINT_MAX(p, rn), mn));

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 14,
                          2 * need * sizeof(mp_limb_t),
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
            fft_small_op_init_borrowed(Y[i] + j, P, fft_small_op_pool_take(pool));
        }
    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&Za, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&Zb, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(T + 0, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(T + 1, P, fft_small_op_pool_take(pool));

    fwd_set_matrices(fj, X, M, rn, Y, M1, mn);
    fwd_set(fj + 8, &A, ap, an);
    fwd_set(fj + 9, &B, bp, bn);
    forward(fj, 10, P, par);

    /* a' = m11 a - m01 b -> Za, b' = m00 b - m10 a -> Zb, and the rows of
       the matrix product with the temporaries T0, T1 */
    buf = fft_small_op_pool_scratch_take(pool, 2 * need, 0);
    out_set_adjust(oj + 0, ap, &Za, Y[1] + 1, &A, Y[0] + 1, &B, nslots_adj, n, p, mn, buf);
    out_set_adjust(oj + 1, bp, &Zb, Y[0] + 0, &B, Y[1] + 0, &A, nslots_adj, n, p, mn, buf + need);
    out_set_matmul(oj + 2, M, X, Y, T + 0, T + 1, nslots_mul, rn, mn);
    outputs(oj, 4, P, par);
    n = adjust_size(n, ap, bp, oj[0].ex[0].hi, oj[1].ex[0].hi);
    matrix_set_size(M, rn, mn);
    fft_small_op_pool_scratch_give(pool, buf);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            fft_small_op_pool_give(pool, X[i][j].data);
            fft_small_op_pool_give(pool, Y[i][j].data);
        }
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_give(pool, Za.data);
    fft_small_op_pool_give(pool, Zb.data);
    fft_small_op_pool_give(pool, T[0].data);
    fft_small_op_pool_give(pool, T[1].data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    return n;
}

/******************************************************************************

    (r; b) <- (a; b) M   (cofactor update in gcdext)

******************************************************************************/

/*
    Sets (r, b) = (m00 a + m10 b, m01 a + m11 b) where a, b have n limbs
    (not necessarily normalized); r and b need room for n + M->n + 1 limbs.
    Returns the new size (normalized), or 0 if fft_small is not used.
*/
mp_size_t
_flint_mpn_hgcd_mul_matrix_vector_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], A, B;
    fwd_job_struct fj[6];
    out_job_struct oj[2];
    mp_size_t mn = M->n, an, bn, zn;
    mp_ptr t;
    ulong nslots;
    int i, j, par;

    if (!use_fft(n, mn, FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF))
        return 0;

    an = bn = n;
    MPN_NORM(ap, an);
    MPN_NORM(bp, bn);
    if (an == 0 || bn == 0)
        return 0;

    R = get_default_mpn_ctx();

    if (!fft_small_plan_init_mpn(P, R, n, mn, 2, 0))
        return 0;

    if (!pool_fits(P, 6))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    nslots = fft_small_plan_mpn_nslots(P, n, mn);
    zn = n + mn + 1;
    par = use_parallel(FLINT_MAX(n, mn));

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 6,
                          2 * zn * sizeof(mp_limb_t),
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));

    fwd_set(fj + 0, &A, ap, an);
    fwd_set(fj + 1, X[0] + 0, M->p[0][0], mn);
    fwd_set(fj + 2, &B, bp, bn);
    fwd_set(fj + 3, X[0] + 1, M->p[0][1], mn);
    fwd_set(fj + 4, X[1] + 0, M->p[1][0], mn);
    fwd_set(fj + 5, X[1] + 1, M->p[1][1], mn);
    forward(fj, 6, P, par);

    /* r = m00 a + m10 b -> X00, b = m01 a + m11 b -> X01. The products
       have at most n + mn + 1 limbs, but like the classical code we may
       only write limb n + mn when it is nonzero (the caller's allocation
       is sized for the normalized result), so the exports go through a
       temporary. */
    t = fft_small_op_pool_scratch_take(pool, 2 * zn, 0);
    for (i = 0; i < 2; i++)
    {
        dot_set2(oj[i].dots.d, X[0] + i, X[0] + i, &A, 0, X[1] + i, &B);
        oj[i].dots.nd = 1;
        exp_set_unsigned(oj[i].ex, t + i * zn, zn, X[0] + i, 2, nslots);
        oj[i].nex = 1;
    }
    outputs(oj, 2, P, par);

    flint_mpn_copyi(rp, t, zn - 1);
    flint_mpn_copyi(bp, t + zn, zn - 1);
    if ((t[zn - 1] | t[2 * zn - 1]) != 0)
    {
        rp[zn - 1] = t[zn - 1];
        bp[zn - 1] = t[2 * zn - 1];
    }
    else
        zn--;

    fft_small_op_pool_scratch_give(pool, t);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_pool_give(pool, X[i][j].data);
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    while ((rp[zn - 1] | bp[zn - 1]) == 0)
        zn--;

    return zn;
}

/******************************************************************************

    (a; b) <- M^{-1} (a; b), result known to fit in nn limbs

******************************************************************************/

/* r = v mod (B^N - 1) in [0, B^N - 1), for v of vn limbs; r may equal v */
static void
fold_bnm1(mp_ptr r, mp_srcptr v, mp_size_t vn, mp_size_t N)
{
    mp_limb_t cy;
    mp_size_t m, i;

    if (vn <= N)
    {
        if (r != v)
            flint_mpn_copyi(r, v, vn);
        flint_mpn_zero(r + vn, N - vn);
    }
    else
    {
        if (r != v)
            flint_mpn_copyi(r, v, N);
        v += N;
        vn -= N;

        while (vn > 0)
        {
            m = FLINT_MIN(vn, N);
            cy = mpn_add(r, r, N, v, m);
            while (cy)
                cy = mpn_add_1(r, r, N, cy);   /* B^N = 1 */
            v += m;
            vn -= m;
        }
    }

    for (i = 0; i < N; i++)
        if (r[i] != UWORD_MAX)
            return;
    flint_mpn_zero(r, N);                      /* B^N - 1 = 0 */
}

int
_flint_mpn_hgcd_matrix_apply_fft(const flint_mpn_hgcd_matrix_struct * M,
    mp_size_t mn[2][2], mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t nn)
{
    mpn_ctx_struct * R;
    fft_small_plan_t P;
    fft_small_op_pool_t pool;
    fft_small_op_struct X[2][2], A, B, T;
    fwd_job_struct fj[7];
    out_job_struct oj[2];
    ulong N, zn, stage;
    mp_ptr t;
    int i, j, par;

    if (!use_fft(nn, M->n, FLINT_MPN_HGCD_APPLY_FFT_CUTOFF))
        return 0;

    R = get_default_mpn_ctx();

    /* The results are nonnegative and less than B^nn, so they are
       determined by their residues mod B^N - 1 for N > nn. Each output
       is a sum of two cyclic products. */
    N = nn + 1;
    if (!fft_small_plan_init_mpn_cyclic(P, R, &N, 2))
        return 0;

    if (!pool_fits(P, 7))
    {
        fft_small_plan_clear(P);
        return 0;
    }

    FLINT_ASSERT(N > (ulong) nn);

    /* the matrix entries must fit (they do in hgcd_reduce), as must a and
       b after one folding */
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            FLINT_ASSERT(mn[i][j] > 0);
    if ((ulong) FLINT_MAX(FLINT_MAX(mn[0][0], mn[0][1]),
                          FLINT_MAX(mn[1][0], mn[1][1])) > N)
    {
        fft_small_plan_clear(P);
        return 0;
    }

    zn = N + (2 * P->bits + P->depth + FLINT_BIT_COUNT(2)) / FLINT_BITS + 1;
    par = use_parallel(nn);

    /* scratch: the folded a and b and B^N - 1 for the forward transforms,
       then N + zn limbs for each of the two exports */
    stage = FLINT_MAX(3 * N, 2 * (N + zn));

    fft_small_op_pool_init(pool);
    fft_small_op_pool_fit(pool, R, fft_small_op_sizeof_data(P), 7,
                          stage * sizeof(mp_limb_t),
                          FFT_SMALL_OP_POOL_RESERVE | FFT_SMALL_OP_POOL_ALLOC);

    t = fft_small_op_pool_scratch_take(pool, stage, 0);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_init_borrowed(X[i] + j, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&A, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&B, P, fft_small_op_pool_take(pool));
    fft_small_op_init_borrowed(&T, P, fft_small_op_pool_take(pool));

    /* a and b mod B^N - 1, and the transform of B^N - 1 */
    fwd_set(fj + 0, &A, ap, n);
    fj[0].fold = 1;
    fj[0].tmp = t;
    fj[0].N = N;
    fwd_set(fj + 1, X[0] + 0, M->p[0][0], mn[0][0]);
    fwd_set(fj + 2, &B, bp, n);
    fj[2].fold = 1;
    fj[2].tmp = t + N;
    fj[2].N = N;
    fwd_set(fj + 3, X[0] + 1, M->p[0][1], mn[0][1]);
    fwd_set(fj + 4, &T, NULL, 0);
    fj[4].ones = 1;
    fj[4].tmp = t + 2 * N;
    fj[4].N = N;
    fwd_set(fj + 5, X[1] + 0, M->p[1][0], mn[1][0]);
    fwd_set(fj + 6, X[1] + 1, M->p[1][1], mn[1][1]);
    forward(fj, 7, P, par);

    /*
        a' = m11 a + m01 ((B^N - 1) - b) = m11 a + m01 T - m01 b -> X11
        b' = m00 b + m10 ((B^N - 1) - a) = m00 b + m10 T - m10 a -> X00
    */
    for (i = 0; i < 2; i++)
    {
        dot_struct * d = oj[i].dots.d;
        exp_job_struct * J = oj[i].ex;
        fft_small_op_struct * Z = (i == 0) ? X[1] + 1 : X[0] + 0;
        fft_small_op_struct * U = (i == 0) ? X[0] + 1 : X[1] + 0;

        d->Z = Z;
        d->X[0] = Z; d->Y[0] = (i == 0) ? &A : &B; d->neg[0] = 0;
        d->X[1] = U; d->Y[1] = &T; d->neg[1] = 0;
        d->X[2] = U; d->Y[2] = (i == 0) ? &B : &A; d->neg[2] = 1;
        d->nterms = 3;
        oj[i].dots.nd = 1;

        J->kind = EXPORT_CYCLIC;
        J->Z = Z;
        J->z = (i == 0) ? ap : bp;
        J->zn = zn;
        J->N = N;
        J->nn = nn;
        J->t = t + i * (N + zn);
        J->buf = J->t + N;
        oj[i].nex = 1;
    }
    outputs(oj, 2, P, par);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            fft_small_op_pool_give(pool, X[i][j].data);
    fft_small_op_pool_give(pool, A.data);
    fft_small_op_pool_give(pool, B.data);
    fft_small_op_pool_give(pool, T.data);
    fft_small_op_pool_scratch_give(pool, t);
    fft_small_op_pool_clear(pool);
    fft_small_plan_clear(P);

    return 1;
}

#endif
