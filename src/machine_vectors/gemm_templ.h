/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
   Internal single-source GEMM implementation, stamped per element type
   (float/double) by gemm.c. All SIMD goes through machine_vectors.h, so
   the same source covers AVX-512 (vec8dz/vec16fz), AVX2 (vec4d/vec8f),
   NEON (vec2d/vec4f), the generic GNU-vector-extension backend and the
   strict ISO C backend with no target #ifdef in this file.

   Expected macros:
     GT_PRE        name prefix (e.g. sg_f32_)
     GT_T          element type (float / double)
     GT_VEC        vector type
     GT_VL         lanes per vector
     GT_MR         microkernel rows (accumulators = GT_MR x 2 vectors)
     GT_LOADU/GT_LOADA/GT_STOREU  vector load/store
     GT_SETZERO/GT_BCAST(ptr)
     GT_FMADD(a,b,c) = a*b + c
     GEMM_PART_SERIAL or GEMM_PART_MT  -- which section to emit

   The serial section emits: _micro, _pack_a, _pack_b, _core.
   The MT section (requires FLINT headers and gemm_mt_barrier from the
   includer) emits: _mt_ctx, _mt_run, _core_mt.

   The includer defines GEMM_PART_SERIAL to emit micro/pack/core, and
   GEMM_PART_MT to emit the thread-pool driver. */

#define GT_CAT_(a, b) a##b
#define GT_CAT(a, b) GT_CAT_(a, b)
#define GT_N(name) GT_CAT(GT_PRE, name)

#define GT_NR (2 * GT_VL)

#ifdef GEMM_PART_SERIAL

/* microkernel: c(GT_MR x GT_NR, row stride ldc) (+)= Ap * Bp over kc */
static void GT_N(micro)(GT_T *c, slong ldc, const GT_T *ap,
                        const GT_T *bp, slong kc, int first)
{
    GT_VEC acc[GT_MR][2];
    for (int r = 0; r < GT_MR; r++)
    {
        acc[r][0] = first ? GT_SETZERO() : GT_LOADU(c + r * ldc);
        acc[r][1] = first ? GT_SETZERO() : GT_LOADU(c + r * ldc + GT_VL);
    }
    for (slong l = 0; l < kc; l++)
    {
        GT_VEC b0 = GT_LOADA(bp + GT_NR * l);
        GT_VEC b1 = GT_LOADA(bp + GT_NR * l + GT_VL);
        for (int r = 0; r < GT_MR; r++)
        {
            GT_VEC a = GT_BCAST(ap + GT_MR * l + r);
            acc[r][0] = GT_FMADD(a, b0, acc[r][0]);
            acc[r][1] = GT_FMADD(a, b1, acc[r][1]);
        }
    }
    for (int r = 0; r < GT_MR; r++)
    {
        GT_STOREU(c + r * ldc, acc[r][0]);
        GT_STOREU(c + r * ldc + GT_VL, acc[r][1]);
    }
}

/* pack (rows x kc) of A into MR-wide l-major panels, zero-filled */
static void GT_N(pack_a)(GT_T *ap, const GT_T *a, slong lda,
                         slong rows, slong kc)
{
    for (slong ir = 0; ir < rows; ir += GT_MR)
    {
        slong rr = rows - ir < GT_MR ? rows - ir : GT_MR;
        for (slong l = 0; l < kc; l++)
        {
            for (slong r = 0; r < rr; r++)
                ap[l * GT_MR + r] = a[(ir + r) * lda + l];
            for (slong r = rr; r < GT_MR; r++)
                ap[l * GT_MR + r] = 0;
        }
        ap += kc * GT_MR;
    }
}

static void GT_N(pack_b)(GT_T *bp, const GT_T *b, slong ldb,
                         slong kc, slong cols)
{
    for (slong jr = 0; jr < cols; jr += GT_NR)
    {
        slong cc = cols - jr < GT_NR ? cols - jr : GT_NR;
        for (slong l = 0; l < kc; l++)
        {
            for (slong j = 0; j < cc; j++)
                bp[l * GT_NR + j] = b[l * ldb + jr + j];
            for (slong j = cc; j < GT_NR; j++)
                bp[l * GT_NR + j] = 0;
        }
        bp += kc * GT_NR;
    }
}

/* thin-shape paths (native element type, non-modular): these shapes
   make panel packing a bad trade -- with m <= MR the B pack dominates
   the row sweep, and with n <= NR the padded microkernel wastes most
   of its width. Both use streaming no-pack kernels instead. */

/* m <= GT_MR: B streamed exactly once, C strips held in registers */
static void GT_N(core_thin_m)(GT_T *C, slong ldc, const GT_T *A,
                              slong lda, const GT_T *B, slong ldb,
                              slong m, slong k, slong n)
{
    slong j0 = 0;
    for (; j0 + 4 * GT_VL <= n; j0 += 4 * GT_VL)
    {
        GT_VEC acc[GT_MR][4];
        for (slong r = 0; r < m; r++)
            for (int v = 0; v < 4; v++) acc[r][v] = GT_SETZERO();
        for (slong l = 0; l < k; l++)
        {
            const GT_T *b = B + l * ldb + j0;
            GT_VEC b0 = GT_LOADU(b), b1 = GT_LOADU(b + GT_VL);
            GT_VEC b2 = GT_LOADU(b + 2 * GT_VL), b3 = GT_LOADU(b + 3 * GT_VL);
            for (slong r = 0; r < m; r++)
            {
                GT_VEC a = GT_BCAST(A + r * lda + l);
                acc[r][0] = GT_FMADD(a, b0, acc[r][0]);
                acc[r][1] = GT_FMADD(a, b1, acc[r][1]);
                acc[r][2] = GT_FMADD(a, b2, acc[r][2]);
                acc[r][3] = GT_FMADD(a, b3, acc[r][3]);
            }
        }
        for (slong r = 0; r < m; r++)
            for (int v = 0; v < 4; v++)
                GT_STOREU(C + r * ldc + j0 + v * GT_VL, acc[r][v]);
    }
    for (; j0 < n; j0++)   /* scalar column tail */
        for (slong r = 0; r < m; r++)
        {
            GT_T sacc = 0;
            for (slong l = 0; l < k; l++)
                sacc += A[r * lda + l] * B[l * ldb + j0];
            C[r * ldc + j0] = sacc;
        }
}

/* n <= GT_NR: B packed once into a tiny k x NR panel; A streamed */
static void GT_N(core_thin_n)(GT_T *C, slong ldc, const GT_T *A,
                              slong lda, const GT_T *B, slong ldb,
                              slong m, slong k, slong n)
{
    GT_T *bp = flint_aligned_alloc(64,
        (k * GT_NR * (slong) sizeof(GT_T) + 63) & ~(slong) 63);
    for (slong l = 0; l < k; l++)
    {
        for (slong j = 0; j < n; j++) bp[l * GT_NR + j] = B[l * ldb + j];
        for (slong j = n; j < GT_NR; j++) bp[l * GT_NR + j] = 0;
    }
    GT_T stage[GT_NR];
    slong i0 = 0;
    for (; i0 + GT_MR <= m; i0 += GT_MR)
    {
        GT_VEC acc[GT_MR][2];
        for (int r = 0; r < GT_MR; r++)
            acc[r][0] = acc[r][1] = GT_SETZERO();
        for (slong l = 0; l < k; l++)
        {
            GT_VEC b0 = GT_LOADA(bp + l * GT_NR);
            GT_VEC b1 = GT_LOADA(bp + l * GT_NR + GT_VL);
            for (int r = 0; r < GT_MR; r++)
            {
                GT_VEC a = GT_BCAST(A + (i0 + r) * lda + l);
                acc[r][0] = GT_FMADD(a, b0, acc[r][0]);
                acc[r][1] = GT_FMADD(a, b1, acc[r][1]);
            }
        }
        for (int r = 0; r < GT_MR; r++)
        {
            GT_STOREU(stage, acc[r][0]);
            GT_STOREU(stage + GT_VL, acc[r][1]);
            memcpy(C + (i0 + r) * ldc, stage, n * sizeof(GT_T));
        }
    }
    for (; i0 < m; i0++)   /* single-row tail */
    {
        GT_VEC a0 = GT_SETZERO(), a1 = GT_SETZERO();
        for (slong l = 0; l < k; l++)
        {
            GT_VEC av = GT_BCAST(A + i0 * lda + l);
            a0 = GT_FMADD(av, GT_LOADA(bp + l * GT_NR), a0);
            a1 = GT_FMADD(av, GT_LOADA(bp + l * GT_NR + GT_VL), a1);
        }
        GT_STOREU(stage, a0);
        GT_STOREU(stage + GT_VL, a1);
        memcpy(C + i0 * ldc, stage, n * sizeof(GT_T));
    }
    flint_aligned_free(bp);
}

#if GEMM_SMALL_DIM > 0

/* all dims <= GEMM_SMALL_DIM (native type, non-modular): direct
   no-pack kernel. The whole problem is L2-resident, so panel packing
   and edge-tile staging only add overhead; instead sweep VL-wide
   column strips with MR-row register accumulators, loading A and B
   in place, with a scalar tail for the last n % VL columns. */
static void GT_N(core_small)(GT_T *C, slong ldc, const GT_T *A,
                             slong lda, const GT_T *B, slong ldb,
                             slong m, slong k, slong n)
{
    slong nv = n / GT_VL * GT_VL;
    slong rem = n - nv;
    /* pad the tail columns once (k <= GEMM_SMALL_DIM, so tiny): the
       vector kernel then covers them too; a scalar tail would be a
       per-element FMA latency chain and dominates at awkward n */
    GT_T btail[GEMM_SMALL_DIM * GT_VL] __attribute__((aligned(64)));
    GT_T stage[GT_VL];
    if (rem)
        for (slong l = 0; l < k; l++)
        {
            for (slong j = 0; j < rem; j++)
                btail[l * GT_VL + j] = B[l * ldb + nv + j];
            for (slong j = rem; j < GT_VL; j++)
                btail[l * GT_VL + j] = 0;
        }
    for (slong i0 = 0; i0 < m; i0 += GT_MR)
    {
        slong rr = m - i0 < GT_MR ? m - i0 : GT_MR;
        for (slong j0 = 0; j0 < nv; j0 += GT_VL)
        {
            if (rr == GT_MR)   /* constant bounds: accs stay in regs */
            {
                GT_VEC acc[GT_MR];
                for (int r = 0; r < GT_MR; r++) acc[r] = GT_SETZERO();
                for (slong l = 0; l < k; l++)
                {
                    GT_VEC b0 = GT_LOADU(B + l * ldb + j0);
                    for (int r = 0; r < GT_MR; r++)
                        acc[r] = GT_FMADD(GT_BCAST(A + (i0 + r) * lda + l),
                                          b0, acc[r]);
                }
                for (int r = 0; r < GT_MR; r++)
                    GT_STOREU(C + (i0 + r) * ldc + j0, acc[r]);
            }
            else               /* single tail row-block per matrix */
            {
                for (slong r = 0; r < rr; r++)
                {
                    GT_VEC a0 = GT_SETZERO();
                    for (slong l = 0; l < k; l++)
                        a0 = GT_FMADD(GT_BCAST(A + (i0 + r) * lda + l),
                                      GT_LOADU(B + l * ldb + j0), a0);
                    GT_STOREU(C + (i0 + r) * ldc + j0, a0);
                }
            }
        }
        if (rem)
        {
            if (rr == GT_MR)
            {
                GT_VEC acc[GT_MR];
                for (int r = 0; r < GT_MR; r++) acc[r] = GT_SETZERO();
                for (slong l = 0; l < k; l++)
                {
                    GT_VEC b0 = GT_LOADA(btail + l * GT_VL);
                    for (int r = 0; r < GT_MR; r++)
                        acc[r] = GT_FMADD(GT_BCAST(A + (i0 + r) * lda + l),
                                          b0, acc[r]);
                }
                for (int r = 0; r < GT_MR; r++)
                {
                    GT_STOREU(stage, acc[r]);
                    memcpy(C + (i0 + r) * ldc + nv, stage,
                           rem * sizeof(GT_T));
                }
            }
            else
                for (slong r = 0; r < rr; r++)
                {
                    GT_VEC a0 = GT_SETZERO();
                    for (slong l = 0; l < k; l++)
                        a0 = GT_FMADD(GT_BCAST(A + (i0 + r) * lda + l),
                                      GT_LOADA(btail + l * GT_VL), a0);
                    GT_STOREU(stage, a0);
                    memcpy(C + (i0 + r) * ldc + nv, stage,
                           rem * sizeof(GT_T));
                }
        }
    }
}

#endif /* GEMM_SMALL_DIM > 0 */

static void GT_N(core)(GT_T *C, slong ldc,
                       const GT_T *A, slong lda,
                       const GT_T *B, slong ldb,
                       slong m, slong k, slong n)
{
    if (k > 0)
    {
        if (m <= GT_MR)
        {
            GT_N(core_thin_m)(C, ldc, A, lda, B, ldb, m, k, n);
            return;
        }
#if GEMM_SMALL_DIM > 0
        if (m <= GEMM_SMALL_DIM && n <= GEMM_SMALL_DIM &&
            k <= GEMM_SMALL_DIM)
        {
            GT_N(core_small)(C, ldc, A, lda, B, ldb, m, k, n);
            return;
        }
#endif
        if (n <= GT_NR)
        {
            GT_N(core_thin_n)(C, ldc, A, lda, B, ldb, m, k, n);
            return;
        }
    }
    slong kcap = k < SG_KC ? k : SG_KC;
    slong ncap = n < SG_NC ? n : SG_NC;
    slong mcap = m < SG_MC ? m : SG_MC;
    slong bpsz = kcap * (ncap + GT_NR) * (slong) sizeof(GT_T);
    slong apsz = kcap * (mcap + GT_MR) * (slong) sizeof(GT_T);
    char *scratch = gemm_get_scratch(bpsz + apsz + 64);
    GT_T *bp = (GT_T *) scratch;
    GT_T *ap = (GT_T *)(scratch + ((bpsz + 63) & ~(slong) 63));
    GT_T stage[GT_MR * GT_NR];

    for (slong jc = 0; jc < n; jc += SG_NC)
    {
        slong nc = n - jc < SG_NC ? n - jc : SG_NC;
        for (slong pc = 0; pc < k; pc += SG_KC)
        {
            slong kc = k - pc < SG_KC ? k - pc : SG_KC;
            int first = pc == 0;
            GT_N(pack_b)(bp, B + pc * ldb + jc, ldb, kc, nc);
            for (slong ic = 0; ic < m; ic += SG_MC)
            {
                slong mc = m - ic < SG_MC ? m - ic : SG_MC;
                GT_N(pack_a)(ap, A + ic * lda + pc, lda, mc, kc);
                for (slong jr = 0; jr < nc; jr += GT_NR)
                {
                    slong cc = nc - jr < GT_NR ? nc - jr : GT_NR;
                    const GT_T *bpp = bp + (jr / GT_NR) * kc * GT_NR;
                    for (slong ir = 0; ir < mc; ir += GT_MR)
                    {
                        slong rr = mc - ir < GT_MR ? mc - ir : GT_MR;
                        const GT_T *app = ap + (ir / GT_MR) * kc * GT_MR;
                        GT_T *ct = C + (ic + ir) * ldc + jc + jr;
                        if (rr == GT_MR && cc == GT_NR)
                            GT_N(micro)(ct, ldc, app, bpp, kc, first);
                        else
                        {
                            if (first)
                                memset(stage, 0, sizeof stage);
                            else
                            {
                                for (slong r = 0; r < rr; r++)
                                    memcpy(stage + r * GT_NR, ct + r * ldc,
                                           cc * sizeof(GT_T));
                                for (slong r = rr; r < GT_MR; r++)
                                    memset(stage + r * GT_NR, 0,
                                           GT_NR * sizeof(GT_T));
                                for (slong r = 0; r < rr; r++)
                                    memset(stage + r * GT_NR + cc, 0,
                                           (GT_NR - cc) * sizeof(GT_T));
                            }
                            GT_N(micro)(stage, GT_NR, app, bpp, kc, 0);
                            for (slong r = 0; r < rr; r++)
                                memcpy(ct + r * ldc, stage + r * GT_NR,
                                       cc * sizeof(GT_T));
                        }
                    }
                }
            }
        }
    }

    gemm_put_scratch(scratch);
}

#endif /* GEMM_PART_SERIAL */

#if defined(GEMM_PART_MT) && defined(GEMM_MT_ATOMIC)

typedef struct
{
    GT_T *C; slong ldc;
    const GT_T *A; slong lda;
    const GT_T *B; slong ldb;
    slong m, k, n;
    int nt;
    slong mc;
    GT_T *bp;
    GT_T **ap;
    _Atomic slong next_jr;
    _Atomic slong next_ic;
    _Atomic int bar_count;
    _Atomic int bar_sense;
} GT_N(mt_ctx);

typedef struct { GT_N(mt_ctx) *ctx; int tid; int sense; } GT_N(mt_arg);

#define GT_BARRIER(g, sensep) \
    gemm_mt_barrier(&(g)->bar_count, &(g)->bar_sense, &(g)->next_jr, \
                    &(g)->next_ic, (g)->nt, sensep)

static void GT_N(mt_run)(void *varg)
{
    GT_N(mt_arg) *w = varg;
    GT_N(mt_ctx) *g = w->ctx;
    GT_T *ap = g->ap[w->tid];
    GT_T stage[GT_MR * GT_NR];

    for (slong jc = 0; jc < g->n; jc += SG_NC)
    {
        slong nc = g->n - jc < SG_NC ? g->n - jc : SG_NC;
        slong npan = (nc + GT_NR - 1) / GT_NR;
        for (slong pc = 0; pc < g->k; pc += SG_KC)
        {
            slong kc = g->k - pc < SG_KC ? g->k - pc : SG_KC;
            int first = pc == 0;

            GT_BARRIER(g, &w->sense);

            for (;;)
            {
                slong jr = atomic_fetch_add(&g->next_jr, 1);
                if (jr >= npan) break;
                slong cols = nc - jr * GT_NR < GT_NR ? nc - jr * GT_NR : GT_NR;
                GT_N(pack_b)(g->bp + jr * kc * GT_NR,
                             g->B + pc * g->ldb + jc + jr * GT_NR,
                             g->ldb, kc, cols);
            }

            GT_BARRIER(g, &w->sense);

            for (;;)
            {
                slong ic = atomic_fetch_add(&g->next_ic, 1) * g->mc;
                if (ic >= g->m) break;
                slong mc = g->m - ic < g->mc ? g->m - ic : g->mc;
                GT_N(pack_a)(ap, g->A + ic * g->lda + pc, g->lda, mc, kc);
                for (slong jr = 0; jr < nc; jr += GT_NR)
                {
                    slong cc = nc - jr < GT_NR ? nc - jr : GT_NR;
                    const GT_T *bpp = g->bp + (jr / GT_NR) * kc * GT_NR;
                    for (slong ir = 0; ir < mc; ir += GT_MR)
                    {
                        slong rr = mc - ir < GT_MR ? mc - ir : GT_MR;
                        const GT_T *app = ap + (ir / GT_MR) * kc * GT_MR;
                        GT_T *ct = g->C + (ic + ir) * g->ldc + jc + jr;
                        if (rr == GT_MR && cc == GT_NR)
                            GT_N(micro)(ct, g->ldc, app, bpp, kc, first);
                        else
                        {
                            if (first)
                                memset(stage, 0, sizeof stage);
                            else
                            {
                                for (slong r = 0; r < rr; r++)
                                    memcpy(stage + r * GT_NR,
                                           ct + r * g->ldc,
                                           cc * sizeof(GT_T));
                                for (slong r = rr; r < GT_MR; r++)
                                    memset(stage + r * GT_NR, 0,
                                           GT_NR * sizeof(GT_T));
                                for (slong r = 0; r < rr; r++)
                                    memset(stage + r * GT_NR + cc, 0,
                                           (GT_NR - cc) * sizeof(GT_T));
                            }
                            GT_N(micro)(stage, GT_NR, app, bpp, kc, 0);
                            for (slong r = 0; r < rr; r++)
                                memcpy(ct + r * g->ldc, stage + r * GT_NR,
                                       cc * sizeof(GT_T));
                        }
                    }
                }
            }

            GT_BARRIER(g, &w->sense);
        }
    }

}

static void GT_N(core_mt)(GT_T *C, slong ldc,
                          const GT_T *A, slong lda,
                          const GT_T *B, slong ldb,
                          slong m, slong k, slong n, slong thread_limit)
{
    /* thin shapes are bandwidth bound; the streaming serial paths win */
    if (m <= GT_MR || n <= GT_NR)
    {
        GT_N(core)(C, ldc, A, lda, B, ldb, m, k, n);
        return;
    }
    if (thread_limit <= 0) thread_limit = flint_get_num_threads();
    double flop = 2.0 * (double) m * (double) n * (double) k;
    slong tcap = (slong)(flop / SG_MT_FLOP_PER_THREAD) + 1;
    if (thread_limit > tcap) thread_limit = tcap;

    thread_pool_handle *handles = NULL;
    slong nw = 0;
    if (thread_limit > 1)
        nw = flint_request_threads(&handles, thread_limit);
    if (nw == 0)
    {
        GT_N(core)(C, ldc, A, lda, B, ldb, m, k, n);
        if (handles != NULL)
            flint_give_back_threads(handles, nw);
        return;
    }
    slong nt = nw + 1;

    GT_N(mt_ctx) g;
    g.C = C; g.ldc = ldc; g.A = A; g.lda = lda; g.B = B; g.ldb = ldb;
    g.m = m; g.k = k; g.n = n;
    g.nt = (int) nt;
    {
        slong tgt = (m + 4 * nt - 1) / (4 * nt);
        tgt = (tgt + GT_MR - 1) / GT_MR * GT_MR;
        if (tgt < GT_MR) tgt = GT_MR;
        if (tgt > SG_MC) tgt = SG_MC;
        g.mc = tgt;
    }
    slong kcap = k < SG_KC ? k : SG_KC;
    slong ncap = n < SG_NC ? n : SG_NC;
    slong mcap = m < SG_MC ? m : SG_MC;
    slong bpsz, apsz;

    /* aligned_alloc requires the size to be a multiple of the alignment */
    bpsz = (kcap * (ncap + GT_NR) * (slong) sizeof(GT_T) + 63) & ~(slong) 63;
    apsz = (kcap * (mcap + GT_MR) * (slong) sizeof(GT_T) + 63) & ~(slong) 63;

    g.bp = flint_aligned_alloc(64, bpsz);
    g.ap = flint_malloc(nt * sizeof(GT_T *));
    for (slong t = 0; t < nt; t++)
        g.ap[t] = flint_aligned_alloc(64, apsz);
    atomic_store(&g.next_jr, 0);
    atomic_store(&g.next_ic, 0);
    atomic_store(&g.bar_count, 0);
    atomic_store(&g.bar_sense, 0);

    GT_N(mt_arg) *args = flint_malloc(nt * sizeof(GT_N(mt_arg)));
    for (slong t = 0; t < nt; t++)
    {
        args[t].ctx = &g;
        args[t].tid = (int) t;
        args[t].sense = 0;
    }

    for (slong t = 0; t < nw; t++)
        thread_pool_wake(global_thread_pool, handles[t], 0,
                         GT_N(mt_run), &args[t]);
    GT_N(mt_run)(&args[nw]);
    for (slong t = 0; t < nw; t++)
        thread_pool_wait(global_thread_pool, handles[t]);

    flint_give_back_threads(handles, nw);
    for (slong t = 0; t < nt; t++)
        flint_aligned_free(g.ap[t]);

    flint_free(g.ap);
    flint_aligned_free(g.bp);
    flint_free(args);
}

#undef GT_BARRIER

#endif /* GEMM_PART_MT */

#if defined(GEMM_PART_MT) && !defined(GEMM_MT_ATOMIC)

/*
    Split-only parallel driver: no atomics and no barriers of any kind.
    The row range of C is divided into one independent sub-problem per
    worker, each of which is an ordinary serial multiplication of a row
    block of A by the whole of B. Blocks are disjoint, so the only
    synchronization is the thread pool's own wake/wait.

    This duplicates the packing of B across workers (an O(k*n) cost per
    worker against O(m*k*n/T) of work per worker), so it is slightly
    slower than the cooperative driver at high thread counts, but it has
    no portability surface beyond the thread pool itself.
*/

typedef struct
{
    GT_T * C; slong ldc;
    const GT_T * A; slong lda;
    const GT_T * B; slong ldb;
    slong m, k, n;
    slong i0, rows;
    int is_worker;
}
GT_N(split_arg);

static void GT_N(split_worker)(void * varg)
{
    GT_N(split_arg) * w = (GT_N(split_arg) *) varg;

    GT_N(core)(w->C + w->i0 * w->ldc, w->ldc,
               w->A + w->i0 * w->lda, w->lda,
               w->B, w->ldb, w->rows, w->k, w->n);

    /*
        Thread pool threads outlive the call, and their packing scratch
        is thread-local, so a worker must release it rather than leave a
        cache attached to a thread that will not run our cleanup
        function. The master keeps its cache for the next call.
    */
    if (w->is_worker)
        gemm_scratch_free_local();
}

static void GT_N(core_split_mt)(GT_T *C, slong ldc,
                                const GT_T *A, slong lda,
                                const GT_T *B, slong ldb,
                                slong m, slong k, slong n,
                                slong thread_limit)
{
    thread_pool_handle * handles;
    GT_N(split_arg) * args;
    slong i, nw, nt, given, pos, tcap;
    double flop;

    /* thin shapes are bandwidth bound; the streaming serial paths win */
    if (m <= GT_MR || n <= GT_NR)
    {
        GT_N(core)(C, ldc, A, lda, B, ldb, m, k, n);
        return;
    }

    /* do not thread work that cannot pay for the hand-off */
    flop = 2.0 * (double) m * (double) n * (double) k;
    tcap = (slong) (flop / SG_MT_FLOP_PER_THREAD) + 1;

    if (thread_limit > tcap)
        thread_limit = tcap;

    if (thread_limit > m)
        thread_limit = m;

    nw = 0;
    handles = NULL;

    if (thread_limit > 1)
        nw = flint_request_threads(&handles, thread_limit);

    if (nw == 0)
    {
        GT_N(core)(C, ldc, A, lda, B, ldb, m, k, n);
        if (handles != NULL)
            flint_give_back_threads(handles, nw);
        return;
    }

    nt = nw + 1;
    args = flint_malloc(nt * sizeof(GT_N(split_arg)));

    /* rows split as evenly as possible; earlier blocks take the
       remainder so the master's block is never the largest */
    pos = 0;

    for (i = 0; i < nt; i++)
    {
        given = m / nt + (i < m % nt ? 1 : 0);

        args[i].C = C; args[i].ldc = ldc;
        args[i].A = A; args[i].lda = lda;
        args[i].B = B; args[i].ldb = ldb;
        args[i].m = m; args[i].k = k; args[i].n = n;
        args[i].i0 = pos;
        args[i].rows = given;
        args[i].is_worker = (i < nw);

        pos += given;
    }

    for (i = 0; i < nw; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0,
                         GT_N(split_worker), &args[i]);

    GT_N(split_worker)(&args[nw]);

    for (i = 0; i < nw; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, nw);
    flint_free(args);
}

#endif /* GEMM_PART_MT */

#undef GT_NR
#undef GT_N
#undef GT_CAT
#undef GT_CAT_
