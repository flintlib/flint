/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Template for a blocked, packed, multithreaded matrix multiplication
    C = A*B over Z/nZ around a register-tile microkernel. Included by
    mul_u32.c, mul_u52.c, mul_fp50.c and mul_k52.c, once per (entry type,
    kernel) pair, after the macros below have been defined; it generates

        BT_NAME(pack_a), BT_NAME(pack_b)   packing into l-major panels
        BT_NAME(micro)                     the MR x NR register tile
        BT_NAME(core)                      blocked serial multiplication
        BT_NAME(core_mt)                   split of C over the thread pool

    (BT_NR is defined here as BT_VL*BT_NACC and undefined at the end; the
    other BT_* parameters are left for the includer to change or undefine).

    Types and geometry

      BT_NAME(x)      prefix for the generated identifiers
      BT_ENTRY        entry type of C, A, B (ulong or uint32_t)
      BT_PACKED       element type of the packed panels
      BT_CTX          modulus context; the generated functions take a
                      const BT_CTX * and only pass it on
      BT_MR           rows of a tile
      BT_VL           lanes of an accumulator vector
      BT_NACC         accumulator vectors per tile row; NR = VL*NACC columns
      BT_KC, BT_MC, BT_NC              blocking sizes, in elements
      BT_MT_MIN_WORK  products (m*k*n) per thread below which the split
                      stops adding threads

    Packing primitives

      BT_LIFT(x, ctx)      entry -> packed element
      BT_BSLOT(j)          slot of column j (0 <= j < NR) in a packed B row

    Microkernel primitives

      BT_ACC, BT_BV, BT_AV     accumulator (VL lanes), B and A operands
      BT_CONSTS                loop-invariant constants derived from ctx
      BT_CONSTS_INIT(ctx)      build them
      BT_ACC_ZERO()
      BT_LOAD_C(p)             VL canonical entries of C -> accumulator
      BT_STORE_C(p, acc)       accumulator holding canonical residues -> C
      BT_LOAD_BSTEP(bv, p)     the NACC B operands of one k step
      BT_LOAD_A(p)             one packed A element as an operand
      BT_MUL_ADD(acc, a, b, C) acc + a*b, lane-wise (C for kernels that
                               reduce each product)
      BT_CADENCE(ctx)          k steps between two BT_FOLD's (>= 1)
      BT_FOLD(acc, C, ctx)     in-loop reduction keeping acc congruent
      BT_FINISH(acc, C, ctx)   accumulator -> canonical residues

    The microkernel computes, for an MR x NR tile c of C with row stride
    ldc, c = Ap*Bp over kc if first is set and c += Ap*Bp otherwise, where
    c holds canonical residues on entry (when not first) and on exit.
*/

#define BT_NR (BT_VL * BT_NACC)

/*
    pack (rows x kc) of A into MR-wide l-major panels, zero-filled:
    ap[l*MR + r] = lift(a[r, l])
*/
static void
BT_NAME(pack_a)(BT_PACKED * ap, const BT_ENTRY * a, slong lda,
                slong rows, slong kc, const BT_CTX * ctx)
{
    slong ir, rr, l, r;

    (void) ctx;     /* not every BT_LIFT looks at it */

    for (ir = 0; ir < rows; ir += BT_MR)
    {
        rr = FLINT_MIN(rows - ir, BT_MR);

        for (l = 0; l < kc; l++)
        {
            for (r = 0; r < rr; r++)
                ap[l * BT_MR + r] = BT_LIFT(a[(ir + r) * lda + l], ctx);
            for (r = rr; r < BT_MR; r++)
                ap[l * BT_MR + r] = 0;
        }

        ap += kc * BT_MR;
    }
}

/*
    pack (kc x cols) of B into NR-wide l-major panels, zero-filled:
    column j of the panel going to slot BT_BSLOT(j) of the row, which is
    the order BT_LOAD_BSTEP wants. Either way accumulator v covers
    columns v*VL .. v*VL + VL - 1 of C, so the tile loads and stores C
    with plain vector accesses.
*/
static void
BT_NAME(pack_b)(BT_PACKED * bp, const BT_ENTRY * b, slong ldb,
                slong kc, slong cols, const BT_CTX * ctx)
{
    slong jr, cc, l, j;

    (void) ctx;     /* not every BT_LIFT looks at it */

    for (jr = 0; jr < cols; jr += BT_NR)
    {
        cc = FLINT_MIN(cols - jr, BT_NR);

        for (l = 0; l < kc; l++)
        {
            const BT_ENTRY * brow = b + l * ldb + jr;
            BT_PACKED * bpl = bp + l * BT_NR;

            for (j = 0; j < BT_NR; j++)
                bpl[j] = 0;
            for (j = 0; j < cc; j++)
                bpl[BT_BSLOT(j)] = BT_LIFT(brow[j], ctx);
        }

        bp += kc * BT_NR;
    }
}

/*
    c (MR x NR, row stride ldc, canonical residues) (+)= Ap * Bp over kc,
    written back as canonical residues
*/
static void
BT_NAME(micro)(BT_ENTRY * c, slong ldc, const BT_PACKED * ap,
               const BT_PACKED * bp, slong kc, int first, const BT_CTX * ctx)
{
    BT_ACC acc[BT_MR][BT_NACC];
    const BT_CONSTS C = BT_CONSTS_INIT(ctx);
    const slong cadence = BT_CADENCE(ctx);
    slong l, stop;
    int r, v;

    for (r = 0; r < BT_MR; r++)
        for (v = 0; v < BT_NACC; v++)
            acc[r][v] = first ? BT_ACC_ZERO()
                              : BT_LOAD_C(c + r * ldc + v * BT_VL);

    l = 0;
    while (l < kc)
    {
        stop = FLINT_MIN(l + cadence, kc);

        /* one iteration is already MR*NACC independent multiply-adds; the
           unrolling of -funroll-loops would only add register pressure */
#if defined(__GNUC__)
# pragma GCC unroll 1
#endif
        for (; l < stop; l++)
        {
            BT_BV bv[BT_NACC];

            BT_LOAD_BSTEP(bv, bp + l * BT_NR);

            for (r = 0; r < BT_MR; r++)
            {
                BT_AV av = BT_LOAD_A(ap + l * BT_MR + r);

                for (v = 0; v < BT_NACC; v++)
                    acc[r][v] = BT_MUL_ADD(acc[r][v], av, bv[v], &C);
            }
        }

        if (l < kc)
            for (r = 0; r < BT_MR; r++)
                for (v = 0; v < BT_NACC; v++)
                    acc[r][v] = BT_FOLD(acc[r][v], &C, ctx);
    }

    for (r = 0; r < BT_MR; r++)
        for (v = 0; v < BT_NACC; v++)
            BT_STORE_C(c + r * ldc + v * BT_VL, BT_FINISH(acc[r][v], &C, ctx));
}

/* blocked serial core ******************************************************/

/* C (m x n, stride ldc) = A (m x k, stride lda) * B (k x n, stride ldb) */
static void
BT_NAME(core)(BT_ENTRY * C, slong ldc, const BT_ENTRY * A, slong lda,
              const BT_ENTRY * B, slong ldb, slong m, slong k, slong n,
              const BT_CTX * ctx)
{
    slong kcap, ncap, mcap, bpsz, apsz;
    slong jc, pc, ic, jr, ir, nc, kc, mc, cc, rr, r;
    char * scratch;
    BT_PACKED * bp, * ap;
    BT_ENTRY stage[BT_MR * BT_NR];

    if (m <= 0 || n <= 0)
        return;

    if (k <= 0)
    {
        for (r = 0; r < m; r++)
            memset(C + r * ldc, 0, n * sizeof(BT_ENTRY));
        return;
    }

    kcap = FLINT_MIN(k, BT_KC);
    ncap = FLINT_MIN(n, BT_NC);
    mcap = FLINT_MIN(m, BT_MC);

    /* aligned_alloc wants sizes that are multiples of the alignment */
    bpsz = (kcap * (ncap + BT_NR) * (slong) sizeof(BT_PACKED) + 63) & ~(slong) 63;
    apsz = (kcap * (mcap + BT_MR) * (slong) sizeof(BT_PACKED) + 63) & ~(slong) 63;
    scratch = flint_aligned_alloc(64, bpsz + apsz);
    bp = (BT_PACKED *) scratch;
    ap = (BT_PACKED *) (scratch + bpsz);

    for (jc = 0; jc < n; jc += BT_NC)
    {
        nc = FLINT_MIN(n - jc, BT_NC);

        for (pc = 0; pc < k; pc += BT_KC)
        {
            int first = (pc == 0);

            kc = FLINT_MIN(k - pc, BT_KC);
            BT_NAME(pack_b)(bp, B + pc * ldb + jc, ldb, kc, nc, ctx);

            for (ic = 0; ic < m; ic += BT_MC)
            {
                mc = FLINT_MIN(m - ic, BT_MC);
                BT_NAME(pack_a)(ap, A + ic * lda + pc, lda, mc, kc, ctx);

                for (jr = 0; jr < nc; jr += BT_NR)
                {
                    const BT_PACKED * bpp = bp + (jr / BT_NR) * kc * BT_NR;

                    cc = FLINT_MIN(nc - jr, BT_NR);

                    for (ir = 0; ir < mc; ir += BT_MR)
                    {
                        const BT_PACKED * app = ap + (ir / BT_MR) * kc * BT_MR;
                        BT_ENTRY * ct = C + (ic + ir) * ldc + jc + jr;

                        rr = FLINT_MIN(mc - ir, BT_MR);

                        if (rr == BT_MR && cc == BT_NR)
                        {
                            BT_NAME(micro)(ct, ldc, app, bpp, kc, first, ctx);
                        }
                        else
                        {
                            /* edge tile: stage through a full tile whose
                               padding rows and columns are zero (and stay
                               zero, since the packed padding is zero) */
                            if (first)
                                memset(stage, 0, sizeof(stage));
                            else
                            {
                                for (r = 0; r < rr; r++)
                                {
                                    memcpy(stage + r * BT_NR, ct + r * ldc,
                                           cc * sizeof(BT_ENTRY));
                                    memset(stage + r * BT_NR + cc, 0,
                                           (BT_NR - cc) * sizeof(BT_ENTRY));
                                }
                                for (r = rr; r < BT_MR; r++)
                                    memset(stage + r * BT_NR, 0,
                                           BT_NR * sizeof(BT_ENTRY));
                            }

                            BT_NAME(micro)(stage, BT_NR, app, bpp, kc, 0, ctx);

                            for (r = 0; r < rr; r++)
                                memcpy(ct + r * ldc, stage + r * BT_NR,
                                       cc * sizeof(BT_ENTRY));
                        }
                    }
                }
            }
        }
    }

    flint_aligned_free(scratch);
}

/* parallel driver **********************************************************/

typedef struct
{
    BT_ENTRY * C; slong ldc;
    const BT_ENTRY * A; slong lda;
    const BT_ENTRY * B; slong ldb;
    slong m, k, n;
    const BT_CTX * ctx;
}
BT_NAME(split_arg);

static void
BT_NAME(split_worker)(void * varg)
{
    BT_NAME(split_arg) * w = (BT_NAME(split_arg) *) varg;
    BT_NAME(core)(w->C, w->ldc, w->A, w->lda, w->B, w->ldb,
                  w->m, w->k, w->n, w->ctx);
}

/*
    Split C into independent blocks along its longer side (rows of A or
    columns of B), one ordinary serial multiplication per thread. Blocks
    are disjoint, so nothing is shared but read-only inputs. The packing
    of the shared operand is duplicated across workers, an O(k*(m+n))
    cost against O(m*k*n/T) per worker.
*/
static void
BT_NAME(core_mt)(BT_ENTRY * C, slong ldc, const BT_ENTRY * A, slong lda,
                 const BT_ENTRY * B, slong ldb, slong m, slong k, slong n,
                 const BT_CTX * ctx, slong thread_limit)
{
    thread_pool_handle * handles = NULL;
    BT_NAME(split_arg) * args;
    slong nw = 0, nt, i, pos, given, len, tcap;
    int split_rows;
    double work;

    work = (double) m * (double) n * (double) k;
    tcap = (slong) (work / BT_MT_MIN_WORK);
    thread_limit = FLINT_MIN(thread_limit, tcap);

    split_rows = (m >= n);
    len = split_rows ? m : n;

    /* each block should be a few tiles */
    thread_limit = FLINT_MIN(thread_limit,
                             len / (split_rows ? 2 * BT_MR : 2 * BT_NR));

    if (thread_limit > 1)
        nw = flint_request_threads(&handles, thread_limit);

    if (nw == 0)
    {
        BT_NAME(core)(C, ldc, A, lda, B, ldb, m, k, n, ctx);
        if (handles != NULL)
            flint_give_back_threads(handles, nw);
        return;
    }

    nt = nw + 1;
    args = flint_malloc(nt * sizeof(BT_NAME(split_arg)));

    pos = 0;
    for (i = 0; i < nt; i++)
    {
        given = len / nt + (i < len % nt ? 1 : 0);

        /* round block boundaries to whole tiles when possible */
        if (i < nt - 1)
        {
            slong tile = split_rows ? BT_MR : BT_NR;
            given = ((given + tile / 2) / tile) * tile;
            given = FLINT_MIN(given, len - pos);
        }
        else
            given = len - pos;

        args[i].ctx = ctx;
        args[i].k = k;

        if (split_rows)
        {
            args[i].C = C + pos * ldc; args[i].ldc = ldc;
            args[i].A = A + pos * lda; args[i].lda = lda;
            args[i].B = B; args[i].ldb = ldb;
            args[i].m = given; args[i].n = n;
        }
        else
        {
            args[i].C = C + pos; args[i].ldc = ldc;
            args[i].A = A; args[i].lda = lda;
            args[i].B = B + pos; args[i].ldb = ldb;
            args[i].m = m; args[i].n = given;
        }

        pos += given;
    }

    for (i = 0; i < nw; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0,
                         BT_NAME(split_worker), &args[i]);

    BT_NAME(split_worker)(&args[nw]);

    for (i = 0; i < nw; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, nw);
    flint_free(args);
}

/* the BT_* parameters are left defined: the includer redefines the ones
   that change between two instantiations */
#undef BT_NR
