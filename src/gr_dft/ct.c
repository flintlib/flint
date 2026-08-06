/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "gr.h"
#include "gr_dft.h"

/* The butterfly loops are organized so that the inner loops perform no
   w-dependent branching: within each block of a pass, the twiddle
   indices j whose roots are special rotations sit at the known
   positions j = 0, hm/4, hm/2, 3 hm/4 (exponents 0, n/8, n/4, 3n/8).
   These few butterflies are handled individually through
   _gr_dft_mul_root, and the remaining ranges of j run through plain
   loops multiplying by a walking pointer into the root table (one
   gr_mul per butterfly), or into the complex Karatsuba tables when
   those are enabled. */

/* DIF butterflies (t = a - b; a += b; b = w^e t) for j in [j0, j1)
   within the block at base index kbase, with generic root
   multiplications. Requires 1 <= j0 <= j1 <= hm and j1 * rstep <= n/2. */
static int
_gr_dft_ct_bflys_dif(gr_ptr x, slong stride, ulong kbase, ulong hm,
        ulong j0, ulong j1, ulong rstep, int inverse,
        gr_ptr t, gr_ptr rtmp, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong xs = stride * sz;
    gr_ptr a = GR_ENTRY(x, (slong) (kbase + j0) * stride, sz);
    gr_ptr b = GR_ENTRY(x, (slong) (kbase + j0 + hm) * stride, sz);
    ulong j, e0 = inverse ? P->n - j0 * rstep : j0 * rstep;
    slong es = (inverse ? -(slong) rstep : (slong) rstep) * sz;
    gr_srcptr w = GR_ENTRY(P->roots, e0, sz);

    if (P->wtab == NULL)
    {
        for (j = j0; j < j1; j++)
        {
            status |= gr_sub(t, a, b, ctx);
            status |= gr_add(a, a, b, ctx);
            status |= gr_mul(b, t, w, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            w = (gr_srcptr) ((char *) w + es);
        }
    }
    else
    {
        gr_ctx_struct * rctx = P->real_ctx;
        slong rsz = rctx->sizeof_elem;
        slong cs = 3 * ((inverse ? -(slong) rstep : (slong) rstep) * rsz);
        gr_srcptr c = GR_ENTRY(P->wtab, 3 * e0, rsz);

        for (j = j0; j < j1; j++)
        {
            gr_srcptr tr = t, ti = GR_ENTRY(t, 1, rsz);
            gr_ptr br = b, bi = GR_ENTRY(b, 1, rsz);

            status |= gr_sub(t, a, b, ctx);
            status |= gr_add(a, a, b, ctx);

            /* 3-multiplication complex Karatsuba product b = w t,
               with (c, d - c, d + c) at the walking table pointer */
            status |= gr_add(rtmp, tr, ti, rctx);
            status |= gr_mul(rtmp, rtmp, c, rctx);
            status |= gr_mul(br, ti, GR_ENTRY(c, 2, rsz), rctx);
            status |= gr_mul(bi, tr, GR_ENTRY(c, 1, rsz), rctx);
            status |= gr_sub(br, rtmp, br, rctx);
            status |= gr_add(bi, rtmp, bi, rctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            c = (gr_srcptr) ((char *) c + cs);
        }
    }

    return status;
}

/* DIT butterflies (t = w^e b; b = a - t; a += t), otherwise as above. */
static int
_gr_dft_ct_bflys_dit(gr_ptr x, slong stride, ulong kbase, ulong hm,
        ulong j0, ulong j1, ulong rstep, int inverse,
        gr_ptr t, gr_ptr rtmp, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong xs = stride * sz;
    gr_ptr a = GR_ENTRY(x, (slong) (kbase + j0) * stride, sz);
    gr_ptr b = GR_ENTRY(x, (slong) (kbase + j0 + hm) * stride, sz);
    ulong j, e0 = inverse ? P->n - j0 * rstep : j0 * rstep;
    slong es = (inverse ? -(slong) rstep : (slong) rstep) * sz;
    gr_srcptr w = GR_ENTRY(P->roots, e0, sz);

    if (P->wtab == NULL)
    {
        for (j = j0; j < j1; j++)
        {
            status |= gr_mul(t, b, w, ctx);
            status |= gr_sub(b, a, t, ctx);
            status |= gr_add(a, a, t, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            w = (gr_srcptr) ((char *) w + es);
        }
    }
    else
    {
        gr_ctx_struct * rctx = P->real_ctx;
        slong rsz = rctx->sizeof_elem;
        slong cs = 3 * ((inverse ? -(slong) rstep : (slong) rstep) * rsz);
        gr_srcptr c = GR_ENTRY(P->wtab, 3 * e0, rsz);

        for (j = j0; j < j1; j++)
        {
            gr_srcptr br = b, bi = GR_ENTRY(b, 1, rsz);
            gr_ptr tr = t, ti = GR_ENTRY(t, 1, rsz);

            /* t = w b by complex Karatsuba */
            status |= gr_add(rtmp, br, bi, rctx);
            status |= gr_mul(rtmp, rtmp, c, rctx);
            status |= gr_mul(tr, bi, GR_ENTRY(c, 2, rsz), rctx);
            status |= gr_mul(ti, br, GR_ENTRY(c, 1, rsz), rctx);
            status |= gr_sub(tr, rtmp, tr, rctx);
            status |= gr_add(ti, rtmp, ti, rctx);

            status |= gr_sub(b, a, t, ctx);
            status |= gr_add(a, a, t, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            c = (gr_srcptr) ((char *) c + cs);
        }
    }

    return status;
}

/* One special butterfly at twiddle index j (a rotation handled through
   the class dispatch of _gr_dft_mul_root). */
static int
_gr_dft_ct_bfly_special(int dit, gr_ptr x, slong stride, ulong kbase,
        ulong hm, ulong j, ulong rstep, int inverse,
        gr_ptr t, gr_ptr rtmp, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_ptr a = GR_ENTRY(x, (slong) (kbase + j) * stride, sz);
    gr_ptr b = GR_ENTRY(x, (slong) (kbase + j + hm) * stride, sz);

    if (dit)
    {
        status |= _gr_dft_mul_root(t, b, j * rstep, inverse, rtmp, P);
        status |= gr_sub(b, a, t, ctx);
        status |= gr_add(a, a, t, ctx);
    }
    else
    {
        status |= gr_sub(t, a, b, ctx);
        status |= gr_add(a, a, b, ctx);
        status |= _gr_dft_mul_root(b, t, j * rstep, inverse, rtmp, P);
    }

    return status;
}

/* Packed-table butterfly loops: the stage-s twiddles are contiguous
   (stage_tab + n - 2^s, entries w^(j rstep), j < hm), so the forward
   walk is sequential; the inverse multiplies by
   w^(-j rstep) = -stage_s[hm - j], a reversed sequential walk with
   the negation folded into the butterfly (swapping the roles of the
   sum and difference), at no extra cost. */
static int
_gr_dft_ct_bflys_dif_packed(gr_ptr x, slong stride, ulong kbase, ulong hm,
        ulong j0, ulong j1, int inverse, gr_srcptr stage,
        gr_ptr t, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong xs = stride * sz;
    gr_ptr a = GR_ENTRY(x, (slong) (kbase + j0) * stride, sz);
    gr_ptr b = GR_ENTRY(x, (slong) (kbase + j0 + hm) * stride, sz);
    ulong j;

    if (!inverse)
    {
        gr_srcptr w = GR_ENTRY(stage, j0, sz);

        for (j = j0; j < j1; j++)
        {
            status |= gr_sub(t, a, b, ctx);
            status |= gr_add(a, a, b, ctx);
            status |= gr_mul(b, t, w, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            w = (gr_srcptr) ((char *) w + sz);
        }
    }
    else
    {
        gr_srcptr w = GR_ENTRY(stage, hm - j0, sz);

        for (j = j0; j < j1; j++)
        {
            status |= gr_sub(t, b, a, ctx);     /* -(a - b) */
            status |= gr_add(a, a, b, ctx);
            status |= gr_mul(b, t, w, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            w = (gr_srcptr) ((char *) w - sz);
        }
    }

    return status;
}

static int
_gr_dft_ct_bflys_dit_packed(gr_ptr x, slong stride, ulong kbase, ulong hm,
        ulong j0, ulong j1, int inverse, gr_srcptr stage,
        gr_ptr t, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong xs = stride * sz;
    gr_ptr a = GR_ENTRY(x, (slong) (kbase + j0) * stride, sz);
    gr_ptr b = GR_ENTRY(x, (slong) (kbase + j0 + hm) * stride, sz);
    ulong j;

    if (!inverse)
    {
        gr_srcptr w = GR_ENTRY(stage, j0, sz);

        for (j = j0; j < j1; j++)
        {
            status |= gr_mul(t, b, w, ctx);
            status |= gr_sub(b, a, t, ctx);
            status |= gr_add(a, a, t, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            w = (gr_srcptr) ((char *) w + sz);
        }
    }
    else
    {
        gr_srcptr w = GR_ENTRY(stage, hm - j0, sz);

        for (j = j0; j < j1; j++)
        {
            status |= gr_mul(t, b, w, ctx);     /* -(b w^(-j rstep)) */
            status |= gr_add(b, a, t, ctx);
            status |= gr_sub(a, a, t, ctx);

            a = (gr_ptr) ((char *) a + xs);
            b = (gr_ptr) ((char *) b + xs);
            w = (gr_srcptr) ((char *) w - sz);
        }
    }

    return status;
}

/* Special rotations at the quarter and eighth positions of a stage
   (complex mode): positional, using only the constant c (the real
   part of w^(n/8), read from the packed entry at hm/4) and sign/swap
   patterns; the inverse rotations are the conjugates. */
static int
_gr_dft_ct_rot_special_packed(gr_ptr r, gr_srcptr v, ulong hm, ulong j,
        int inverse, gr_srcptr stage, gr_ptr rtmp,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * rctx = P->real_ctx;
    slong rsz = rctx->sizeof_elem;
    gr_srcptr xr = v, xi = GR_ENTRY(v, 1, rsz);
    gr_ptr rr = r, ri = GR_ENTRY(r, 1, rsz);

    if (2 * j == hm || hm == 2)
    {
        /* w^(n/4) = -i (forward), w^(-n/4) = i (inverse) */
        if (!inverse)
        {
            status |= gr_set(rr, xi, rctx);
            status |= gr_neg(ri, xr, rctx);
        }
        else
        {
            status |= gr_neg(rr, xi, rctx);
            status |= gr_set(ri, xr, rctx);
        }
    }
    else
    {
        gr_srcptr c = GR_ENTRY(stage, hm / 4, ctx->sizeof_elem);

        if (4 * j == hm)
        {
            if (!inverse)
            {
                /* c (1 - i): rr = c (xr + xi), ri = c (xi - xr) */
                status |= gr_add(rtmp, xr, xi, rctx);
                status |= gr_mul(rr, c, rtmp, rctx);
                status |= gr_sub(rtmp, xi, xr, rctx);
                status |= gr_mul(ri, c, rtmp, rctx);
            }
            else
            {
                /* c (1 + i): rr = c (xr - xi), ri = c (xr + xi) */
                status |= gr_sub(rtmp, xr, xi, rctx);
                status |= gr_mul(rr, c, rtmp, rctx);
                status |= gr_add(rtmp, xr, xi, rctx);
                status |= gr_mul(ri, c, rtmp, rctx);
            }
        }
        else    /* j = 3 hm / 4 */
        {
            if (!inverse)
            {
                /* -c (1 + i): rr = c (xi - xr), ri = -c (xr + xi) */
                status |= gr_sub(rtmp, xi, xr, rctx);
                status |= gr_mul(rr, c, rtmp, rctx);
                status |= gr_add(rtmp, xr, xi, rctx);
                status |= gr_mul(ri, c, rtmp, rctx);
                status |= gr_neg(ri, ri, rctx);
            }
            else
            {
                /* -c (1 - i): rr = -c (xr + xi), ri = c (xr - xi) */
                status |= gr_add(rtmp, xr, xi, rctx);
                status |= gr_mul(rr, c, rtmp, rctx);
                status |= gr_neg(rr, rr, rctx);
                status |= gr_sub(rtmp, xr, xi, rctx);
                status |= gr_mul(ri, c, rtmp, rctx);
            }
        }
    }

    return status;
}

static int
_gr_dft_ct_bfly_special_packed(int dit, gr_ptr x, slong stride, ulong kbase,
        ulong hm, ulong j, int inverse, gr_srcptr stage,
        gr_ptr t, gr_ptr rtmp, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_ptr a = GR_ENTRY(x, (slong) (kbase + j) * stride, sz);
    gr_ptr b = GR_ENTRY(x, (slong) (kbase + j + hm) * stride, sz);

    if (dit)
    {
        status |= _gr_dft_ct_rot_special_packed(t, b, hm, j, inverse,
                stage, rtmp, P, ctx);
        status |= gr_sub(b, a, t, ctx);
        status |= gr_add(a, a, t, ctx);
    }
    else
    {
        gr_ptr u = t;   /* difference, then rotated in place is unsafe;
                           rotate from the freshly computed difference
                           held in t into b directly */
        status |= gr_sub(t, a, b, ctx);
        status |= gr_add(a, a, b, ctx);
        status |= _gr_dft_ct_rot_special_packed(b, u, hm, j, inverse,
                stage, rtmp, P, ctx);
    }

    return status;
}

/* The butterflies of one block with twiddle indices in [jlo, jhi):
   the free rotation at j = 0 (when jlo = 0), the special rotations at
   j = hm/4, hm/2, 3 hm/4 (in complex mode), and plain loops in
   between. The full block is [0, hm); sub-ranges are used by the
   threaded driver. */
static int
_gr_dft_ct_jrange(int dit, gr_ptr x, slong stride, ulong kbase, ulong hm,
        ulong jlo, ulong jhi, ulong rstep, int inverse,
        gr_ptr t, gr_ptr rtmp, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong cur;

    /* j = 0: (a, b) -> (a + b, a - b) in both directions */
    if (jlo == 0)
    {
        gr_ptr a = GR_ENTRY(x, (slong) kbase * stride, ctx->sizeof_elem);
        gr_ptr b = GR_ENTRY(x, (slong) (kbase + hm) * stride, ctx->sizeof_elem);

        status |= gr_sub(t, a, b, ctx);
        status |= gr_add(a, a, b, ctx);
        status |= gr_set(b, t, ctx);
    }

    cur = FLINT_MAX(jlo, 1);
    if (cur >= jhi)
        return status;

    if (P->stage_tab != NULL)
    {
        gr_srcptr stage = GR_ENTRY(P->stage_tab,
                P->n - 2 * hm, ctx->sizeof_elem);

        if (P->real_ctx == NULL)
        {
            if (dit)
                status |= _gr_dft_ct_bflys_dit_packed(x, stride, kbase, hm,
                        cur, jhi, inverse, stage, t, P, ctx);
            else
                status |= _gr_dft_ct_bflys_dif_packed(x, stride, kbase, hm,
                        cur, jhi, inverse, stage, t, P, ctx);
        }
        else
        {
            ulong sp[3];
            slong nsp = 0, i;

            if (hm == 2)
            {
                sp[nsp++] = 1;
            }
            else if (hm >= 4)
            {
                ulong q = hm / 4;
                sp[nsp++] = q;
                sp[nsp++] = 2 * q;
                sp[nsp++] = 3 * q;
            }

            for (i = 0; i < nsp && cur < jhi; i++)
            {
                ulong sj = sp[i];

                if (sj < cur)
                    continue;

                if (cur < FLINT_MIN(sj, jhi))
                {
                    if (dit)
                        status |= _gr_dft_ct_bflys_dit_packed(x, stride,
                                kbase, hm, cur, FLINT_MIN(sj, jhi),
                                inverse, stage, t, P, ctx);
                    else
                        status |= _gr_dft_ct_bflys_dif_packed(x, stride,
                                kbase, hm, cur, FLINT_MIN(sj, jhi),
                                inverse, stage, t, P, ctx);
                }

                if (sj < jhi)
                    status |= _gr_dft_ct_bfly_special_packed(dit, x, stride,
                            kbase, hm, sj, inverse, stage, t, rtmp, P, ctx);

                cur = FLINT_MAX(cur, sj + 1);
            }

            if (cur < jhi)
            {
                if (dit)
                    status |= _gr_dft_ct_bflys_dit_packed(x, stride, kbase,
                            hm, cur, jhi, inverse, stage, t, P, ctx);
                else
                    status |= _gr_dft_ct_bflys_dif_packed(x, stride, kbase,
                            hm, cur, jhi, inverse, stage, t, P, ctx);
            }
        }

        return status;
    }

    if (P->wclass == NULL)
    {
        if (dit)
            status |= _gr_dft_ct_bflys_dit(x, stride, kbase, hm, cur, jhi,
                    rstep, inverse, t, rtmp, P, ctx);
        else
            status |= _gr_dft_ct_bflys_dif(x, stride, kbase, hm, cur, jhi,
                    rstep, inverse, t, rtmp, P, ctx);
    }
    else
    {
        /* special twiddle positions within the block */
        ulong sp[3];
        slong nsp = 0, i;

        if (hm == 2)
        {
            sp[nsp++] = 1;
        }
        else if (hm >= 4)
        {
            ulong q = hm / 4;
            sp[nsp++] = q;
            sp[nsp++] = 2 * q;
            sp[nsp++] = 3 * q;
        }

        for (i = 0; i < nsp && cur < jhi; i++)
        {
            ulong sj = sp[i];

            if (sj < cur)
                continue;

            if (cur < FLINT_MIN(sj, jhi))
            {
                if (dit)
                    status |= _gr_dft_ct_bflys_dit(x, stride, kbase, hm,
                            cur, FLINT_MIN(sj, jhi), rstep, inverse,
                            t, rtmp, P, ctx);
                else
                    status |= _gr_dft_ct_bflys_dif(x, stride, kbase, hm,
                            cur, FLINT_MIN(sj, jhi), rstep, inverse,
                            t, rtmp, P, ctx);
            }

            if (sj < jhi)
                status |= _gr_dft_ct_bfly_special(dit, x, stride, kbase,
                        hm, sj, rstep, inverse, t, rtmp, P, ctx);

            cur = FLINT_MAX(cur, sj + 1);
        }

        if (cur < jhi)
        {
            if (dit)
                status |= _gr_dft_ct_bflys_dit(x, stride, kbase, hm,
                        cur, jhi, rstep, inverse, t, rtmp, P, ctx);
            else
                status |= _gr_dft_ct_bflys_dif(x, stride, kbase, hm,
                        cur, jhi, rstep, inverse, t, rtmp, P, ctx);
        }
    }

    return status;
}

/* Decimation in frequency: natural order input, bit-reversed output. */
/* All decimation levels s, ..., 1 over the aligned block
   [kbase, kbase + 2^s), in a depth-first blocked order: blocks larger
   than GR_DFT_CT_BLOCK_BYTES run their top pass, then recurse into
   the two halves, so that each cache-sized sub-block completes all
   its remaining levels while resident; the level-by-level order,
   which would stream the block from memory once per level, is used
   only once a block fits. The set of butterflies and their
   dependency order are identical to the plain level-by-level
   transform, so results are bitwise the same. */
static int
_gr_dft_ct_dif_range(gr_ptr x, slong stride, ulong kbase, int s,
        int inverse, gr_ptr t, gr_ptr rtmp,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong n = P->n;

    if (s >= 2 &&
        (UWORD(1) << s) * (ulong) ctx->sizeof_elem > GR_DFT_CT_BLOCK_BYTES)
    {
        ulong hm = UWORD(1) << (s - 1);

        status |= _gr_dft_ct_jrange(0, x, stride, kbase, hm, 0, hm,
                n >> s, inverse, t, rtmp, P, ctx);
        status |= _gr_dft_ct_dif_range(x, stride, kbase, s - 1,
                inverse, t, rtmp, P, ctx);
        status |= _gr_dft_ct_dif_range(x, stride, kbase + hm, s - 1,
                inverse, t, rtmp, P, ctx);
    }
    else
    {
        int ss;
        ulong k, m;

        for (ss = s; ss >= 1; ss--)
        {
            m = UWORD(1) << ss;
            for (k = kbase; k < kbase + (UWORD(1) << s); k += m)
                status |= _gr_dft_ct_jrange(0, x, stride, k, m >> 1,
                        0, m >> 1, n >> ss, inverse, t, rtmp, P, ctx);
        }
    }

    return status;
}

static int
_gr_dft_ct_dit_range(gr_ptr x, slong stride, ulong kbase, int s,
        int inverse, gr_ptr t, gr_ptr rtmp,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong n = P->n;

    if (s >= 2 &&
        (UWORD(1) << s) * (ulong) ctx->sizeof_elem > GR_DFT_CT_BLOCK_BYTES)
    {
        ulong hm = UWORD(1) << (s - 1);

        status |= _gr_dft_ct_dit_range(x, stride, kbase, s - 1,
                inverse, t, rtmp, P, ctx);
        status |= _gr_dft_ct_dit_range(x, stride, kbase + hm, s - 1,
                inverse, t, rtmp, P, ctx);
        status |= _gr_dft_ct_jrange(1, x, stride, kbase, hm, 0, hm,
                n >> s, inverse, t, rtmp, P, ctx);
    }
    else
    {
        int ss;
        ulong k, m;

        for (ss = 1; ss <= s; ss++)
        {
            m = UWORD(1) << ss;
            for (k = kbase; k < kbase + (UWORD(1) << s); k += m)
                status |= _gr_dft_ct_jrange(1, x, stride, k, m >> 1,
                        0, m >> 1, n >> ss, inverse, t, rtmp, P, ctx);
        }
    }

    return status;
}

static int
_gr_dft_ct_dif(gr_ptr x, slong stride, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t, rtmp = NULL;

    GR_TMP_INIT(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    status |= _gr_dft_ct_dif_range(x, stride, 0, P->depth,
            inverse, t, rtmp, P, ctx);

    GR_TMP_CLEAR(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    return status;
}

/* Decimation in time: bit-reversed input, natural order output. */
static int
_gr_dft_ct_dit(gr_ptr x, slong stride, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t, rtmp = NULL;

    GR_TMP_INIT(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    status |= _gr_dft_ct_dit_range(x, stride, 0, P->depth,
            inverse, t, rtmp, P, ctx);

    GR_TMP_CLEAR(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    return status;
}

/* In-place transform of the strided vector x[0], x[stride], ...,
   x[(n-1)*stride]. If scrambled is set, the forward transform leaves
   the output in bit-reversed order and the inverse transform expects
   its input in bit-reversed order; otherwise both are in natural order.
   The inverse transform omits the scaling by 1/n. */
int
_gr_dft_ct(gr_ptr x, slong stride, int inverse, int scrambled,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (P->depth == 0)
        return GR_SUCCESS;

    if (!inverse)
    {
        status |= _gr_dft_ct_dif(x, stride, 0, P, ctx);
        if (!scrambled)
            _gr_dft_bit_reverse(x, stride, P->depth, ctx);
    }
    else
    {
        if (!scrambled)
            _gr_dft_bit_reverse(x, stride, P->depth, ctx);
        status |= _gr_dft_ct_dit(x, stride, 1, P, ctx);
    }

    return status;
}

/* Pass-parallel threaded execution (used for top-level transforms
   only; the serial _gr_dft_ct above is also invoked for the
   sub-transforms inside the Bailey workers and must never attempt to
   acquire threads itself).

   Within each pass all n/2 butterflies are independent; they are
   indexed by flat units u in [0, n/2), u -> (block u / hm, twiddle
   index u mod hm), so that contiguous unit ranges map to runs of
   whole blocks with clipped ranges at the ends. This distributes the
   work evenly both in the early passes of a DIF transform (one block
   of n/2 butterflies, split by twiddle ranges) and in the late ones
   (many small blocks). Each pass is one scatter/join round, giving
   log2(n) synchronization points, in exchange for full parallel
   width at any size (in contrast with the Bailey algorithm, whose
   three phases have only n1 or n2 work items). */

typedef struct
{
    gr_ptr x;
    slong stride;
    int dit;
    int inverse;
    ulong m;
    ulong rstep;
    ulong lo;
    ulong hi;
    const gr_dft_pre_struct * P;
    gr_ctx_struct * ctx;
    int status;
}
_gr_dft_ct_work_t;

static void
_gr_dft_ct_worker(void * arg)
{
    _gr_dft_ct_work_t * w = arg;
    gr_ctx_struct * ctx = w->ctx;
    const gr_dft_pre_struct * P = w->P;
    int status = GR_SUCCESS;
    ulong hm = w->m >> 1;
    ulong kb, kb0 = w->lo / hm, kb1 = (w->hi - 1) / hm;
    gr_ptr t, rtmp = NULL;

    GR_TMP_INIT(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    for (kb = kb0; kb <= kb1; kb++)
    {
        ulong jlo = (kb == kb0) ? w->lo - kb0 * hm : 0;
        ulong jhi = (kb == kb1) ? w->hi - kb1 * hm : hm;

        status |= _gr_dft_ct_jrange(w->dit, w->x, w->stride, kb * w->m,
                hm, jlo, jhi, w->rstep, w->inverse, t, rtmp, P, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    w->status = status;
}

static int
_gr_dft_ct_pass_mt(int dit, gr_ptr x, slong stride, int s, int inverse,
        slong nchunks, thread_pool_handle * handles,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong n = P->n, units = n / 2;
    slong i;
    _gr_dft_ct_work_t * args;

    args = flint_malloc(nchunks * sizeof(_gr_dft_ct_work_t));

    for (i = 0; i < nchunks; i++)
    {
        args[i].x = x;
        args[i].stride = stride;
        args[i].dit = dit;
        args[i].inverse = inverse;
        args[i].m = UWORD(1) << s;
        args[i].rstep = n >> s;
        args[i].lo = (units * (ulong) i) / (ulong) nchunks;
        args[i].hi = (units * (ulong) (i + 1)) / (ulong) nchunks;
        args[i].P = P;
        args[i].ctx = ctx;
        args[i].status = GR_SUCCESS;

        if (i < nchunks - 1)
            thread_pool_wake(global_thread_pool, handles[i], 0,
                    _gr_dft_ct_worker, &args[i]);
    }

    _gr_dft_ct_worker(&args[nchunks - 1]);

    for (i = 0; i < nchunks - 1; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < nchunks; i++)
        status |= args[i].status;

    flint_free(args);
    return status;
}

/* Owner-computes tail: after the top t decimation levels, the array
   decomposes into 2^t independent contiguous sub-blocks of size
   2^(depth - t); a worker runs ALL the remaining passes over its own
   range of blocks without further synchronization, keeping the data
   in its private cache. (The pass-parallel top levels are the only
   ones whose butterflies genuinely span the array; sweeping every
   level across all workers instead is limited by the shared memory
   bandwidth, which does not grow with the thread count.) */
typedef struct
{
    gr_ptr x;
    slong stride;
    int dit;
    int inverse;
    int sdepth;             /* depth of one sub-block */
    ulong blo;
    ulong bhi;              /* range of sub-block indices */
    const gr_dft_pre_struct * P;
    gr_ctx_struct * ctx;
    int status;
}
_gr_dft_ct_tail_work_t;

static void
_gr_dft_ct_tail_worker(void * arg)
{
    _gr_dft_ct_tail_work_t * w = arg;
    gr_ctx_struct * ctx = w->ctx;
    const gr_dft_pre_struct * P = w->P;
    int status = GR_SUCCESS;
    ulong b;
    ulong bm = UWORD(1) << w->sdepth;
    gr_ptr t, rtmp = NULL;

    GR_TMP_INIT(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    if (!w->dit)
    {
        for (b = w->blo; b < w->bhi; b++)
            status |= _gr_dft_ct_dif_range(w->x, w->stride, b * bm,
                    w->sdepth, w->inverse, t, rtmp, P, ctx);
    }
    else
    {
        for (b = w->blo; b < w->bhi; b++)
            status |= _gr_dft_ct_dit_range(w->x, w->stride, b * bm,
                    w->sdepth, w->inverse, t, rtmp, P, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    w->status = status;
}

static int
_gr_dft_ct_tail_mt(int dit, gr_ptr x, slong stride, int sdepth, int inverse,
        slong nchunks, thread_pool_handle * handles,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    ulong nblocks = P->n >> sdepth;
    slong i;
    _gr_dft_ct_tail_work_t * args;

    args = flint_malloc(nchunks * sizeof(_gr_dft_ct_tail_work_t));

    for (i = 0; i < nchunks; i++)
    {
        args[i].x = x;
        args[i].stride = stride;
        args[i].dit = dit;
        args[i].inverse = inverse;
        args[i].sdepth = sdepth;
        args[i].blo = (nblocks * (ulong) i) / (ulong) nchunks;
        args[i].bhi = (nblocks * (ulong) (i + 1)) / (ulong) nchunks;
        args[i].P = P;
        args[i].ctx = ctx;
        args[i].status = GR_SUCCESS;

        if (i < nchunks - 1)
            thread_pool_wake(global_thread_pool, handles[i], 0,
                    _gr_dft_ct_tail_worker, &args[i]);
    }

    _gr_dft_ct_tail_worker(&args[nchunks - 1]);

    for (i = 0; i < nchunks - 1; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < nchunks; i++)
        status |= args[i].status;

    flint_free(args);
    return status;
}

/* Top-level Cooley-Tukey transform with pass-parallel threading when
   worker threads are available (falling back to the serial transform
   otherwise). Must only be used for top-level transforms. */
int
_gr_dft_ct_threaded(gr_ptr x, slong stride, int inverse, int scrambled,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong blk, max_chunks, nchunks;
    int threadsafe, s;
    thread_pool_handle * handles = P->threads;
    slong num_workers = P->num_threads;
    thread_pool_handle * pool_handles = NULL;
    slong pool_workers = 0;

    if (P->depth == 0)
        return GR_SUCCESS;

    threadsafe = (gr_ctx_is_threadsafe(ctx) == T_TRUE) &&
            (P->real_ctx == NULL ||
             gr_ctx_is_threadsafe(P->real_ctx) == T_TRUE);

    if (!threadsafe)
    {
        handles = NULL;
        num_workers = 0;
    }

    /* threading granularity: at most n / serial_block chunks per pass */
    blk = (P->serial_block > 0) ? P->serial_block : GR_DFT_SERIAL_BLOCK_DEFAULT;
    max_chunks = (slong) (P->n / (ulong) blk);

    if (threadsafe && handles == NULL && max_chunks >= 2 &&
        flint_get_num_threads() > 1)
    {
        pool_workers = flint_request_threads(&pool_handles,
                FLINT_MIN(max_chunks, (slong) (P->n / 2)));
        handles = pool_handles;
        num_workers = pool_workers;
    }

    if (handles == NULL || num_workers < 1 || max_chunks < 2)
    {
        if (pool_handles != NULL)
            flint_give_back_threads(pool_handles, pool_workers);
        return _gr_dft_ct(x, stride, inverse, scrambled, P, ctx);
    }

    nchunks = FLINT_MIN(num_workers + 1, max_chunks);
    nchunks = FLINT_MIN(nchunks, (slong) (P->n / 2));
    nchunks = FLINT_MAX(nchunks, 1);

    /* the top t levels have butterflies spanning the array and run
       pass-parallel; the remaining levels decompose into 2^t
       independent blocks handled owner-computes by the workers */
    {
        int t = 0;
        while ((slong) (UWORD(1) << t) < nchunks && t < P->depth)
            t++;

        if (!inverse)
        {
            for (s = P->depth; s >= P->depth - t + 1; s--)
                status |= _gr_dft_ct_pass_mt(0, x, stride, s, 0,
                        nchunks, handles, P, ctx);
            if (P->depth - t >= 1)
                status |= _gr_dft_ct_tail_mt(0, x, stride, P->depth - t, 0,
                        nchunks, handles, P, ctx);
            if (!scrambled)
                _gr_dft_bit_reverse(x, stride, P->depth, ctx);
        }
        else
        {
            if (!scrambled)
                _gr_dft_bit_reverse(x, stride, P->depth, ctx);
            if (P->depth - t >= 1)
                status |= _gr_dft_ct_tail_mt(1, x, stride, P->depth - t, 1,
                        nchunks, handles, P, ctx);
            for (s = P->depth - t + 1; s <= P->depth; s++)
                status |= _gr_dft_ct_pass_mt(1, x, stride, s, 1,
                        nchunks, handles, P, ctx);
        }
    }

    if (pool_handles != NULL)
        flint_give_back_threads(pool_handles, pool_workers);

    return status;
}
