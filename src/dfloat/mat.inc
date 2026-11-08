/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* ---- matrix multiplication (plain types) ----

   Included from template.inc once per format, after the vector code.

   For N = 1 the product is flint_dgemm.  For N > 1 it is a SIMD kernel
   over blocks of DF_MW (4, or 8 with AVX-512) columns of the result:
   the columns of B are packed by blocks into structure-of-arrays form
   (component c of the entries (kk, j0 .. j0 + DF_MW - 1) is one
   vector), and each row of A then runs along a block, each entry of
   the row broadcast to the lanes.  The product of a row entry and a
   packed entry is not renormalized: its order sums T[0 .. N - 1] (as
   in the multiplication kernel: exact TwoProds and TwoSums except for
   the plain additions of the last order, relative to the term) go into
   accumulators S[0 .. N] of the order sums of the whole dot product,
   each addition into S[0 .. N - 1] a TwoSum whose residual goes on to
   the next order, the additions into S[N] plain.  The accumulators are
   renormalized every DFLOAT_MM_RENORM terms (which keeps them ordered,
   so that the residuals stay at their order) and at the end.  Dropped,
   per term, as in the multiplication: the low parts of the order N - 1
   products, the products of order >= N and the rounding errors of the
   last order of T; per addition into S[N], its rounding error, which
   is of order N + 1 relative to the sum.  Without the extra order
   S[N], the additions into S[N - 1] would round relative to the
   running sum rather than to the term (losing 5 bits for d4 at k = 64
   on random data); with it, the result is more accurate than a dot
   product of the renormalized products, at about two thirds of the
   operations.

   The complex product is the real product of the m x 2k matrix of the
   parts of A by the 2k x 2n matrix [Re B, Im B; -Im B, Re B] (with rows
   and columns interleaved as the parts are), whose result is the m x 2n
   matrix of the parts of C: four real terms per complex term, with no
   copy of A or C. */

#if defined(DFLOAT_HAVE_SIMD) && N > 1

/* The kernel's vectors: 8 lanes with AVX-512, else the vec4d of the
   rest of the module.  Only exact operations (with a fused fmsub for
   the TwoProd) are used, so the results do not depend on the width. */
#if defined(__AVX512F__) && !defined(FLINT_MACHINE_VECTORS_GENERIC) && !defined(DFLOAT_MM_NO_AVX512)
#define DF_MV vec8dz
#define DF_MW 8
#define DF_MV_ADD(a, b) vec8dz_add(a, b)
#define DF_MV_SUB(a, b) vec8dz_sub(a, b)
#define DF_MV_MUL(a, b) vec8dz_mul(a, b)
#define DF_MV_FMS(a, b, c) vec8dz_fmsub(a, b, c)   /* always fused */
#define DF_MV_SET_D(x) vec8dz_set_d(x)
#define DF_MV_ZERO() vec8dz_zero()
#define DF_MV_LOAD(p) vec8dz_load_unaligned(p)
#define DF_MV_STORE(p, a) vec8dz_store_unaligned(p, a)
#else
#define DF_MV vec4d
#define DF_MW 4
#define DF_MV_ADD(a, b) vec4d_add(a, b)
#define DF_MV_SUB(a, b) vec4d_sub(a, b)
#define DF_MV_MUL(a, b) vec4d_mul(a, b)
#define DF_MV_FMS(a, b, c) vec4d_fmsub_fused(a, b, c)
#define DF_MV_SET_D(x) vec4d_set_d(x)
#define DF_MV_ZERO() vec4d_zero()
#define DF_MV_LOAD(p) vec4d_load_unaligned(p)
#define DF_MV_STORE(p, a) vec4d_store_unaligned(p, a)
#endif

#ifndef DFLOAT_MM_RENORM
#define DFLOAT_MM_RENORM 32
#endif

#ifndef DFLOAT_MM_ROWS
#define DFLOAT_MM_ROWS 2
#endif

FLINT_FORCE_INLINE void
DF(_mm_two_sum)(DF_MV * s, DF_MV * e, DF_MV a, DF_MV b)
{
    DF_MV ss = DF_MV_ADD(a, b);
    DF_MV bb = DF_MV_SUB(ss, a);
    *e = DF_MV_ADD(DF_MV_SUB(a, DF_MV_SUB(ss, bb)), DF_MV_SUB(b, bb));
    *s = ss;
}

/* The accumulators: the order sums S[0 .. N] (one more than the
   number of components, so that the additions into S[N - 1] are exact
   too and only those into S[N] round). */
#define DF_MM_S (N + 1)

/* add t to the order sum o, the residuals rippling down the orders */
FLINT_FORCE_INLINE void
DF(_mm_add)(DF_MV * S, DF_MV t, int o)
{
    int j;
    UNROLL
    for (j = o; j < DF_MM_S - 1; j++)
        DF(_mm_two_sum)(&S[j], &t, S[j], t);
    S[DF_MM_S - 1] = DF_MV_ADD(S[DF_MM_S - 1], t);
}

/* S += a b, a broadcast from the N components at a, b in the lanes:
   the order sums T of the product, as in the multiplication kernel
   (exact except for the plain additions of the last order, whose
   terms are of order N - 1 relative to the product), each then added
   to the accumulator of its order */
FLINT_FORCE_INLINE void
DF(_mm_term)(DF_MV * S, const double * a, const DF_MV * b)
{
    DF_MV av[N], h[N][N], l[N][N], T[N], acc;
    DF_MV R0[N * N], R1[N * N];
    DF_MV * rprev = R0, * rnext = R1, * tmp;
    int i, j, k, nprev = 0, nnext;

    UNROLL
    for (i = 0; i < N; i++)
        av[i] = DF_MV_SET_D(a[i]);

    UNROLL
    for (i = 0; i < N; i++)
    {
        UNROLL
        for (j = 0; i + j < N; j++)
        {
            h[i][j] = DF_MV_MUL(av[i], b[j]);
            if (i + j < N - 1)
                l[i][j] = DF_MV_FMS(av[i], b[j], h[i][j]);
        }
    }

    T[0] = h[0][0];
    UNROLL
    for (k = 1; k < N; k++)
    {
        nnext = 0;
        acc = h[0][k];
#define DF_MM_ACC(x) \
        do { \
            if (k == N - 1) \
                acc = DF_MV_ADD(acc, (x)); \
            else \
                DF(_mm_two_sum)(&acc, &rnext[nnext++], acc, (x)); \
        } while (0)
        UNROLL
        for (i = 1; i <= k; i++)
            DF_MM_ACC(h[i][k - i]);
        UNROLL
        for (i = 0; i < k; i++)
            DF_MM_ACC(l[i][k - 1 - i]);
        UNROLL
        for (i = 0; i < nprev; i++)
            DF_MM_ACC(rprev[i]);
#undef DF_MM_ACC
        T[k] = acc;
        tmp = rprev; rprev = rnext; rnext = tmp;
        nprev = nnext;
    }

    UNROLL
    for (k = 0; k < N; k++)
        DF(_mm_add)(S, T[k], k);
}

/* VecSum and VecSumErr over the order sums, as in the kernels */
FLINT_FORCE_INLINE void
DF(_mm_renorm)(DF_MV * S)
{
    int i;
    UNROLL
    for (i = DF_MM_S - 2; i >= 0; i--)
        DF(_mm_two_sum)(&S[i], &S[i + 1], S[i], S[i + 1]);
    UNROLL
    for (i = 0; i < DF_MM_S - 1; i++)
        DF(_mm_two_sum)(&S[i], &S[i + 1], S[i], S[i + 1]);
}

/* entry (kk, j) of B, or for cplx of [Re B, Im B; -Im B, Re B] with
   rows and columns interleaved (B then has complex entries, ldb in
   complex entries) */
FLINT_FORCE_INLINE void
DF(_mm_bget)(double * r, DF(srcptr) B, slong ldb, slong kk, slong j, int cplx)
{
    int c;
    if (!cplx)
    {
        for (c = 0; c < N; c++)
            r[c] = B[kk * ldb + j].d[c];
    }
    else
    {
        /* the parts of B_{kk/2, j/2}: entries 2 e (re) and 2 e + 1 (im) */
        slong e = (kk / 2) * ldb + j / 2;
        int s = kk & 1, t = j & 1;
        DF(srcptr) p = B + 2 * e + ((s == t) ? 0 : 1);
        int neg = (s == 1 && t == 0);
        for (c = 0; c < N; c++)
            r[c] = neg ? -p->d[c] : p->d[c];
    }
}

/* the k x DF_MW block of columns j0 .. j0 + nc - 1 of B (zeros
   beyond), in SoA form: Bp[(kk N + c) DF_MW + lane] */
static void
DF(_mm_pack)(double * Bp, DF(srcptr) B, slong ldb, slong k, slong j0, slong nc, int cplx)
{
    slong kk, l;
    int c;
    double t[N];

    if (!cplx && nc == DF_MW)
    {
        for (kk = 0; kk < k; kk++)
        {
            DF(srcptr) row = B + kk * ldb + j0;
            double * out = Bp + kk * N * DF_MW;
            for (l = 0; l < DF_MW; l++)
                for (c = 0; c < N; c++)
                    out[DF_MW * c + l] = row[l].d[c];
        }
        return;
    }

    for (kk = 0; kk < k; kk++)
    {
        double * out = Bp + kk * N * DF_MW;
        for (l = 0; l < DF_MW; l++)
        {
            if (l < nc)
                DF(_mm_bget)(t, B, ldb, kk, j0 + l, cplx);
            else
                for (c = 0; c < N; c++)
                    t[c] = 0.0;
            for (c = 0; c < N; c++)
                out[DF_MW * c + l] = t[c];
        }
    }
}

/* entry (i, j) of the product in scalar, for the nonfinite results */
static void
DF(_mm_entry_scalar)(DF(ptr) res, DF(srcptr) A, slong lda, DF(srcptr) B, slong ldb, slong k, slong i, slong j, int cplx)
{
    DF(t) s, b, t;
    slong kk;
    DF(_zero_inline)(s);
    for (kk = 0; kk < k; kk++)
    {
        DF(_mm_bget)(b->d, B, ldb, kk, j, cplx);
        DF(_mul_impl)(t, A + i * lda + kk, b);
        DF(_add_impl)(s, s, t);
    }
    DF(_set_inline)(res, s);
}

/* rows i .. i + R - 1 of the block of columns j0 .. j0 + nc - 1 */
FLINT_FORCE_INLINE void
DF(_mm_rows)(DF(ptr) C, slong ldc, DF(srcptr) A, slong lda, const double * Bp,
    slong k, slong i, slong j0, slong nc, int R,
    DF(srcptr) B, slong ldb, int cplx)
{
    DF_MV S[DFLOAT_MM_ROWS][DF_MM_S];
    slong kk, kend;
    int r, c, l, q;

    UNROLL
    for (r = 0; r < R; r++)
    {
        UNROLL
        for (c = 0; c < DF_MM_S; c++)
            S[r][c] = DF_MV_ZERO();
    }

    for (kk = 0; kk < k; kk = kend)
    {
        kend = FLINT_MIN(k, kk + DFLOAT_MM_RENORM);
        for ( ; kk < kend; kk++)
        {
            DF_MV b[N];
            const double * bp = Bp + kk * N * DF_MW;
            UNROLL
            for (c = 0; c < N; c++)
                b[c] = DF_MV_LOAD(bp + DF_MW * c);
            UNROLL
            for (r = 0; r < R; r++)
                DF(_mm_term)(S[r], A[(i + r) * lda + kk].d, b);
        }
        UNROLL
        for (r = 0; r < R; r++)
            DF(_mm_renorm)(S[r]);
    }

    for (r = 0; r < R; r++)
    {
        double t[N][DF_MW];

        /* (S[N] is below the last component after the renormalization) */
        UNROLL
        for (c = 0; c < N; c++)
            DF_MV_STORE(t[c], S[r][c]);

        for (q = 0; q < DF_MW && q < nc; q += 4)
        {
            vec4d v[N], chk;
            int bad;

            UNROLL
            for (c = 0; c < N; c++)
                v[c] = vec4d_load_unaligned(t[c] + q);
            DF(_vshift_zero_head)(v);
            chk = DF(_vabs_sum)(v);
            bad = vec4d_movemask(vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)));
            UNROLL
            for (c = 0; c < N; c++)
                vec4d_store_unaligned(t[c] + q, v[c]);

            for (l = q; l < q + 4 && l < nc; l++)
            {
                DF(ptr) out = C + (i + r) * ldc + j0 + l;
                if (FLINT_UNLIKELY((bad >> (l - q)) & 1))
                    DF(_mm_entry_scalar)(out, A, lda, B, ldb, k, i + r, j0 + l, cplx);
                else
                {
                    UNROLL
                    for (c = 0; c < N; c++)
                        out->d[c] = t[c][l];
                }
            }
        }
    }
}

/* C (m x n) = A (m x k) B (k x n), row strides in elements; C must not
   overlap A or B.  For cplx, A and C are the matrices of parts (m x k
   and m x n with k and n even, strides in real elements) and B is the
   complex k/2 x n/2 matrix (stride in complex entries).

   B is packed once, block by block of DF_MW columns.  The rows of A
   are taken by panels of mc rows small enough to stay in the L2 cache
   (about 256 KB) while every block of columns runs along them; the
   panels are independent tasks, run in parallel with the threads
   FLINT is set to use. */

typedef struct
{
    DF(ptr) C;
    DF(srcptr) A;
    DF(srcptr) B;
    const double * Bp;
    slong ldc, lda, ldb, m, k, n, mc;
    int cplx;
}
DF(_mm_args_struct);

static void
DF(_mm_worker)(slong t, void * vargs)
{
    DF(_mm_args_struct) * args = vargs;
    slong i, i0, i1, j0, nc, k = args->k;

    i0 = t * args->mc;
    i1 = FLINT_MIN(args->m, i0 + args->mc);

    for (j0 = 0; j0 < args->n; j0 += DF_MW)
    {
        const double * Bp = args->Bp + (j0 / DF_MW) * k * N * DF_MW;
        nc = FLINT_MIN(DF_MW, args->n - j0);

        i = i0;
#if DFLOAT_MM_ROWS > 1
        for ( ; i + DFLOAT_MM_ROWS <= i1; i += DFLOAT_MM_ROWS)
            DF(_mm_rows)(args->C, args->ldc, args->A, args->lda, Bp, k, i, j0, nc, DFLOAT_MM_ROWS, args->B, args->ldb, args->cplx);
#endif
        for ( ; i < i1; i++)
            DF(_mm_rows)(args->C, args->ldc, args->A, args->lda, Bp, k, i, j0, nc, 1, args->B, args->ldb, args->cplx);
    }
}

static void
DFV(_mat_mul_simd)(DF(ptr) C, slong ldc, DF(srcptr) A, slong lda,
    DF(srcptr) B, slong ldb, slong m, slong k, slong n, int cplx)
{
    DF(_mm_args_struct) args;
    double * Bp;
    slong j0, nblocks, ntasks, mc;

    nblocks = (n + DF_MW - 1) / DF_MW;
    Bp = flint_malloc(sizeof(double) * DF_MW * N * FLINT_MAX(k, 1) * nblocks);

    for (j0 = 0; j0 < n; j0 += DF_MW)
        DF(_mm_pack)(Bp + (j0 / DF_MW) * k * N * DF_MW, B, ldb, k, j0, FLINT_MIN(DF_MW, n - j0), cplx);

    /* rows per panel: about 256 KB of A, a multiple of the rows of the
       kernel, and at least two panels per thread when threaded */
    mc = (32768 / (N * FLINT_MAX(k, 1))) / DFLOAT_MM_ROWS * DFLOAT_MM_ROWS;
    mc = FLINT_MAX(mc, DFLOAT_MM_ROWS);

    args.C = C; args.A = A; args.B = B; args.Bp = Bp;
    args.ldc = ldc; args.lda = lda; args.ldb = ldb;
    args.m = m; args.k = k; args.n = n; args.cplx = cplx;

    if (flint_get_num_threads() > 1 && (double) m * k * n * N * N > 1e6)
    {
        slong nt = flint_get_num_threads();
        while (mc > DFLOAT_MM_ROWS && (m + mc - 1) / mc < 2 * nt)
            mc = FLINT_MAX(DFLOAT_MM_ROWS, (mc / 2) / DFLOAT_MM_ROWS * DFLOAT_MM_ROWS);
        args.mc = mc;
        ntasks = (m + mc - 1) / mc;
        flint_parallel_do(DF(_mm_worker), &args, ntasks, -1, FLINT_PARALLEL_STRIDED);
    }
    else
    {
        slong t;
        args.mc = mc;
        ntasks = (m + mc - 1) / mc;
        for (t = 0; t < ntasks; t++)
            DF(_mm_worker)(t, &args);
    }

    flint_free(Bp);
}

#undef DF_MV
#undef DF_MW
#undef DF_MV_ADD
#undef DF_MV_SUB
#undef DF_MV_MUL
#undef DF_MV_FMS
#undef DF_MV_SET_D
#undef DF_MV_ZERO
#undef DF_MV_LOAD
#undef DF_MV_STORE
#undef DF_MM_S

#endif

/* nonzero if the n doubles from x and the m doubles from y overlap */
FLINT_FORCE_INLINE int
DF(_mm_overlap)(const void * x, slong n, const void * y, slong m)
{
    const char * a = x, * b = y;
    return a < b + m * sizeof(double) && b < a + n * sizeof(double);
}

/* the span in doubles of the r x c matrix with stride ld (elements of
   w doubles) */
FLINT_FORCE_INLINE slong
DF(_mm_span)(slong r, slong c, slong ld, slong w)
{
    return (r == 0 || c == 0) ? 0 : ((r - 1) * ld + c) * w;
}

/* Above these dimensions the product is arb_mat_approx_mul or
   acb_mat_approx_mul at 53 N + 4 bits on the exact values (the block
   algorithms of arb, through exact integer matrix products, become
   faster than the kernel, whose cost is cubic; the cutoffs include the
   conversions to and from arf, a quarter of the time at the cutoff for
   d4; with complex entries that are all real, acb is faster still,
   which is not taken into account).  Their result is accurate
   in the norm of the blocks rather than entrywise.

   This makes dfloat depend on arb_mat and acb_mat, which it otherwise
   sits below; it is a provisional solution until there is a fast
   arbitrary-precision floating-point matrix multiplication at a lower
   level. (The approximate functions of arb_mat and acb_mat call dfloat
   for LU decompositions and solving, but never for approx_mul, so there
   is no cycle.) */
#ifndef DFLOAT_MAT_MUL_ARB_CUTOFF
#define DFLOAT_MAT_MUL_ARB_CUTOFF(cplx) \
    ((cplx) ? ((N == 2) ? 1500 : (N == 3) ? 900 : 500) \
            : ((N == 2) ? 1600 : (N == 3) ? 1000 : 600))
#endif

#if N > 1
static int
DF(_mm_all_finite)(const gr_mat_t A, int w)
{
    slong i, j;
    for (i = 0; i < A->r; i++)
    {
        const double * row = ((const double *) A->entries) + i * A->stride * w;
        for (j = 0; j < A->c * w; j += N)
            if (!(fabs(row[j]) <= DBL_MAX))
                return 0;
    }
    return 1;
}

/* (not static: the test calls it directly, below the cutoffs) */
void DFV(_mat_mul_arb)(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cplx);

void
DFV(_mat_mul_arb)(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cplx)
{
    slong i, j, prec = 53 * N + 4;
    slong lda = A->stride, ldb = B->stride, ldc = C->stride;
    const double * Ad = A->entries, * Bd = B->entries;
    double * Cd = C->entries;

    if (!cplx)
    {
        arb_mat_t RA, RB, RC;
        arb_mat_init(RA, A->r, A->c);
        arb_mat_init(RB, B->r, B->c);
        arb_mat_init(RC, C->r, C->c);
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                _dfloat_get_arf(arb_midref(arb_mat_entry(RA, i, j)), Ad + (i * lda + j) * N, N);
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                _dfloat_get_arf(arb_midref(arb_mat_entry(RB, i, j)), Bd + (i * ldb + j) * N, N);
        arb_mat_approx_mul(RC, RA, RB, prec);
        for (i = 0; i < C->r; i++)
            for (j = 0; j < C->c; j++)
                _dfloat_set_arf(Cd + (i * ldc + j) * N, N, NULL, arb_midref(arb_mat_entry(RC, i, j)));
        arb_mat_clear(RA);
        arb_mat_clear(RB);
        arb_mat_clear(RC);
    }
    else
    {
        acb_mat_t RA, RB, RC;
        acb_mat_init(RA, A->r, A->c);
        acb_mat_init(RB, B->r, B->c);
        acb_mat_init(RC, C->r, C->c);
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
            {
                const double * x = Ad + (i * lda + j) * 2 * N;
                _dfloat_get_arf(arb_midref(acb_realref(acb_mat_entry(RA, i, j))), x, N);
                _dfloat_get_arf(arb_midref(acb_imagref(acb_mat_entry(RA, i, j))), x + N, N);
            }
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
            {
                const double * x = Bd + (i * ldb + j) * 2 * N;
                _dfloat_get_arf(arb_midref(acb_realref(acb_mat_entry(RB, i, j))), x, N);
                _dfloat_get_arf(arb_midref(acb_imagref(acb_mat_entry(RB, i, j))), x + N, N);
            }
        acb_mat_approx_mul(RC, RA, RB, prec);
        for (i = 0; i < C->r; i++)
            for (j = 0; j < C->c; j++)
            {
                double * x = Cd + (i * ldc + j) * 2 * N;
                _dfloat_set_arf(x, N, NULL, arb_midref(acb_realref(acb_mat_entry(RC, i, j))));
                _dfloat_set_arf(x + N, N, NULL, arb_midref(acb_imagref(acb_mat_entry(RC, i, j))));
            }
        acb_mat_clear(RA);
        acb_mat_clear(RB);
        acb_mat_clear(RC);
    }
}
#endif

/* The product for the gr methods: C = A B with A and B real (cplx = 0)
   or complex (cplx = 1) matrices of dN parts.  Returns GR_UNABLE where
   there is no fast code (no SIMD, N > 1), for the caller to use the
   classical product. */
static int
DFV(_mat_mul_gr)(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, int cplx, gr_ctx_t ctx)
{
    slong m = A->r, k = A->c, n = B->c, w = cplx ? 2 * N : N;
    slong lda = A->stride, ldb = B->stride, ldc = C->stride;
    const double * Ad = (const double *) A->entries, * Bd = (const double *) B->entries;
    double * Cd = (double *) C->entries;

    if (k != B->r || m != C->r || n != C->c)
        return GR_DOMAIN;

    if (m == 0 || n == 0)
        return GR_SUCCESS;

    if (k == 0)
        return gr_mat_zero(C, ctx);

    /* the result through a temporary if it overlaps an operand */
    if (DF(_mm_overlap)(Cd, DF(_mm_span)(m, n, ldc, w), Ad, DF(_mm_span)(m, k, lda, w))
        || DF(_mm_overlap)(Cd, DF(_mm_span)(m, n, ldc, w), Bd, DF(_mm_span)(k, n, ldb, w)))
    {
        gr_mat_t T;
        int status;
        gr_mat_init(T, m, n, ctx);
        status = DFV(_mat_mul_gr)(T, A, B, cplx, ctx);
        if (status == GR_SUCCESS)
            status = gr_mat_set(C, T, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

#if N == 1
    if (!cplx)
    {
        flint_dgemm(m, k, n, Ad, lda, Bd, ldb, Cd, ldc);
    }
    else
    {
        /* [Re B, Im B; -Im B, Re B], interleaved, as a 2k x 2n matrix */
        double * Bt = flint_malloc(sizeof(double) * 4 * k * n);
        slong kk, j;
        for (kk = 0; kk < k; kk++)
        {
            const double * row = Bd + 2 * kk * ldb;
            double * r0 = Bt + (2 * kk) * (2 * n);
            double * r1 = r0 + 2 * n;
            for (j = 0; j < n; j++)
            {
                double re = row[2 * j], im = row[2 * j + 1];
                r0[2 * j] = re;
                r0[2 * j + 1] = im;
                r1[2 * j] = -im;
                r1[2 * j + 1] = re;
            }
        }
        flint_dgemm(m, 2 * k, 2 * n, Ad, 2 * lda, Bt, 2 * n, Cd, 2 * ldc);
        flint_free(Bt);
    }
    return GR_SUCCESS;
#else
    if (FLINT_MIN(m, FLINT_MIN(k, n)) >= DFLOAT_MAT_MUL_ARB_CUTOFF(cplx)
        && DF(_mm_all_finite)(A, w) && DF(_mm_all_finite)(B, w))
    {
        DFV(_mat_mul_arb)(C, A, B, cplx);
        return GR_SUCCESS;
    }
#if defined(DFLOAT_HAVE_SIMD)
    if (!cplx)
        DFV(_mat_mul_simd)((DF(ptr)) Cd, ldc, (DF(srcptr)) Ad, lda, (DF(srcptr)) Bd, ldb, m, k, n, 0);
    else
        DFV(_mat_mul_simd)((DF(ptr)) Cd, 2 * ldc, (DF(srcptr)) Ad, 2 * lda, (DF(srcptr)) Bd, ldb, m, 2 * k, 2 * n, 1);
    return GR_SUCCESS;
#else
    return GR_UNABLE;
#endif
#endif
}

/* ---- LU decomposition and triangular solving: the basecases ----

   The block recursive algorithms of gr_mat do most of the work in
   gr_mat_mul; what they leave to the basecases is the factorization of
   tall panels of a few columns and the solution of small triangular
   systems with many right-hand sides.  The generic classical
   algorithms work on rows, which are short there, and pay for a
   function call or a scalar operation per entry; the basecases here
   are organized for the shapes they get.

   The LU decomposition is the classical one (the same pivoting and the
   same layout of L as gr_mat_lu_classical, with the pivot of largest
   head) with the elimination done column by column, four rows at a
   time in the SIMD lanes: e_j = a_{j,col} / pivot for the rows below
   the pivot, then column by column x_{j,c} -= e_j p_c with p_c broadcast.
   The results are those of the scalar operations (the same products
   and differences; a block with a nonfinite result is redone in
   scalar). */

#define DF_ENTRY(M, i, j) (((DF(ptr)) (M)->entries) + (i) * (M)->stride + (j))

/* the row in start_row .. m - 1 with the largest head in the column,
   -1 if the column is zero there (a nan head is taken at once) */
static slong
DF(_lu_pivot)(const gr_mat_t M, slong start_row, slong col)
{
    slong i, best = -1;
    double bv = 0.0, v;
    for (i = start_row; i < M->r; i++)
    {
        v = DF_ENTRY(M, i, col)->d[0];
        if (v != v)
            return i;
        v = fabs(v);
        if (v > bv)
        {
            bv = v;
            best = i;
        }
    }
    return best;
}

/* rows j0 .. m - 1 below the pivot row (row, col): the multipliers into
   column rank_col, the remaining columns updated */
static void
DF(_lu_eliminate)(gr_mat_t LU, slong row, slong col, slong rank_col, DF(srcptr) d)
{
    slong m = LU->r, n = LU->c, j = row + 1, c;
    DF(struct) e[4];
    DF(t) t;
    DF(srcptr) prow = DF_ENTRY(LU, row, 0);

#ifdef DFLOAT_HAVE_SIMD
    int stride = (int) (LU->stride * N);
    for ( ; j + 4 <= m; j += 4)
    {
        slong k;
        vec4d ev[N];
        for (k = 0; k < 4; k++)
            DF(_mul_impl)(e + k, DF_ENTRY(LU, j + k, col), d);
        DF(_vload4)(ev, e[0].d, N);
        for (c = col + 1; c < n; c++)
        {
            vec4d pv[N], xv[N], tv[N], rv[N];
            DF(ptr) x = DF_ENTRY(LU, j, c);
            int i;
            for (i = 0; i < N; i++)
                pv[i] = vec4d_set_d(prow[c].d[i]);
            DF(_vload4)(xv, x->d, stride);
            DF(_v_mul_kernel)(tv, NULL, ev, pv, 0, DFLOAT_KERNEL_STRONG);
            DF(_vshift_zero_head)(tv);
            DF(_v_addsub_kernel)(rv, NULL, xv, tv, 1, 0, DFLOAT_KERNEL_STRONG);
            DF(_vshift_zero_head)(rv);
            if (FLINT_LIKELY(!VANY(vec4d_cmp_nle(vec4d_add(DF(_vabs_sum)(tv), DF(_vabs_sum)(rv)), vec4d_set_d(DBL_MAX)))))
            {
                DF(_vstore4)(x->d, rv, stride);
            }
            else
            {
                for (k = 0; k < 4; k++)
                {
                    DF(_mul_impl)(t, e + k, prow + c);
                    DF(_sub_impl)(DF_ENTRY(LU, j + k, c), DF_ENTRY(LU, j + k, c), t);
                }
            }
        }
        for (k = 0; k < 4; k++)
        {
            DF(_zero_inline)(DF_ENTRY(LU, j + k, col));
            DF(_set_inline)(DF_ENTRY(LU, j + k, rank_col), e + k);
        }
    }
#endif

    for ( ; j < m; j++)
    {
        DF(_mul_impl)(e, DF_ENTRY(LU, j, col), d);
        for (c = col + 1; c < n; c++)
        {
            DF(_mul_impl)(t, e, prow + c);
            DF(_sub_impl)(DF_ENTRY(LU, j, c), DF_ENTRY(LU, j, c), t);
        }
        DF(_zero_inline)(DF_ENTRY(LU, j, col));
        DF(_set_inline)(DF_ENTRY(LU, j, rank_col), e);
    }
}

int DFV(gr_mat_lu_classical)(slong * res_rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);

int
DFV(gr_mat_lu_classical)(slong * res_rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong i, m = A->r, n = A->c, r, rank = 0, row = 0, col = 0;
    DF(t) d;
    int status = GR_SUCCESS;

    for (i = 0; i < m; i++)
        P[i] = i;

    if (m == 0 || n == 0)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    status |= gr_mat_set(LU, A, ctx);

    while (row < m && col < n)
    {
        r = DF(_lu_pivot)(LU, row, col);

        if (r == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }
            col++;
            continue;
        }

        rank++;

        if (r != row)
            _gr_mat_swap_rows(LU, P, row, r, ctx);

        DF(inv)(d, DF_ENTRY(LU, row, col));
        DF(_lu_eliminate)(LU, row, col, rank - 1, d);

        row++;
        col++;
    }

    *res_rank = rank;
    return status;
}

/* Triangular solving with many right-hand sides: the rows of X as
   vectors, X_j = (B_j - sum_i L_ji X_i) / L_jj by vector products by a
   scalar (in the order i = 0, 1, ... for tril, n - 1, n - 2, ... for
   triu), with the SIMD vector operations of the ring; for few
   right-hand sides, the classical algorithm of gr_mat (dot products
   along the rows of L). */
static int
DF(_mat_solve_tri_rows)(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, int upper, gr_ctx_t ctx)
{
    slong n = L->r, m = B->c, i, j, jj, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr inv;

    if (X != B)
        status |= gr_mat_set(X, B, ctx);

    GR_TMP_INIT(inv, ctx);

    for (jj = 0; jj < n; jj++)
    {
        j = upper ? n - 1 - jj : jj;
        gr_ptr Xj = GR_MAT_ENTRY(X, j, 0, sz);

        if (upper)
        {
            for (i = j + 1; i < n; i++)
                status |= _gr_vec_submul_scalar(Xj, GR_MAT_ENTRY(X, i, 0, sz), m, GR_MAT_ENTRY(L, j, i, sz), ctx);
        }
        else
        {
            for (i = 0; i < j; i++)
                status |= _gr_vec_submul_scalar(Xj, GR_MAT_ENTRY(X, i, 0, sz), m, GR_MAT_ENTRY(L, j, i, sz), ctx);
        }

        if (!unit)
        {
            status |= gr_inv(inv, GR_MAT_ENTRY(L, j, j, sz), ctx);
            status |= _gr_vec_mul_scalar(Xj, Xj, m, inv, ctx);
        }

        if (status != GR_SUCCESS)
            break;
    }

    GR_TMP_CLEAR(inv, ctx);
    return status;
}

#undef DF_ENTRY
