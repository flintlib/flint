/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* Mixed-radix decimation in time: with len = p m and p prime,

       X_(q + m r) = sum_b v^(b r) (w_len^(b q) Y_b(q)),   v = w_len^m,

   where Y_b is the transform of length m of the b-th residue class
   of the input. The recursion transforms the p residue classes into
   the segments res[b m, (b+1) m), then for each column q combines the
   twiddled values z_b with a transform of prime length p.

   The prime transform uses the principality identity
   1 + v + ... + v^(p-1) = 0 to eliminate the exponent p - 1: writing
   e_b = b r mod p and b* for the index with e_(b*) = p - 1,

       Y_r = sum_b v^(e_b) z_b = sum_(b != b*) v^(e_b) (z_b - z_(b*)),

   which requires p - 2 root multiplications per output r != 0 (the
   b = 0 term has e_0 = 0), hence (p-1)(p-2) in total instead of the
   naive (p-1)^2. If a Bluestein sub-plan for the prime is present
   (P->P1), it is used instead of the direct kernel.

   All root exponents are taken in the table of the top-level plan P:
   at a level with sub-length len, the root of the sub-transform is
   w^rstep with rstep = n / len, twiddles use exponents b q rstep < n,
   and the prime kernel uses exponents e_b (m rstep) = e_b (n / p) < n,
   so no reductions are needed. */

static int
_mixed(gr_ptr res, gr_srcptr vec, slong stride, ulong rstep, ulong len,
        slong lvl, int inverse, gr_ptr ztmp, gr_ptr ytmp,
        gr_ptr t1, gr_ptr rtmp, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    ulong p, m, b, q, r, e, bstar, kstep;

    if (len == 1)
        return gr_set(res, vec, ctx);

    p = P->radices[lvl];
    m = len / p;
    kstep = m * rstep;   /* the transform of length p uses the root w^kstep */

    for (b = 0; b < p; b++)
        status |= _mixed(GR_ENTRY(res, b * m, sz), GR_ENTRY(vec, b * stride, sz),
                p * stride, p * rstep, m, lvl + 1, inverse,
                ztmp, ytmp, t1, rtmp, P, ctx);

    for (q = 0; q < m; q++)
    {
        /* load and twiddle: z_b = w^(b q rstep) res[b m + q] */
        status |= gr_set(ztmp, GR_ENTRY(res, q, sz), ctx);
        for (b = 1, e = 0; b < p; b++)
        {
            e += q * rstep;
            status |= _gr_dft_mul_root(GR_ENTRY(ztmp, b, sz),
                    GR_ENTRY(res, b * m + q, sz), e, inverse, rtmp, P);
        }

        if (P->P1 != NULL && P->P1->n == p)
        {
            /* transform of length p by the Bluestein sub-plan */
            status |= _gr_dft_bluestein(ytmp, ztmp, inverse, P->P1, ctx);
            for (r = 0; r < p; r++)
                status |= gr_set(GR_ENTRY(res, r * m + q, sz),
                        GR_ENTRY(ytmp, r, sz), ctx);
        }
        else
        {
            /* direct kernel; the slots res[r m + q] have all been read */
            gr_ptr acc = GR_ENTRY(res, q, sz);

            status |= gr_set(acc, ztmp, ctx);
            for (b = 1; b < p; b++)
                status |= gr_add(acc, acc, GR_ENTRY(ztmp, b, sz), ctx);

            for (r = 1; r < p; r++)
            {
                /* find b* with b* r = p - 1 mod p */
                for (b = 1, e = r; e != p - 1; b++)
                {
                    e += r;
                    if (e >= p)
                        e -= p;
                }
                bstar = b;

                acc = GR_ENTRY(res, r * m + q, sz);
                status |= gr_sub(acc, ztmp, GR_ENTRY(ztmp, bstar, sz), ctx);

                for (b = 1, e = 0; b < p; b++)
                {
                    e += r;
                    if (e >= p)
                        e -= p;

                    if (b == bstar)
                        continue;

                    status |= gr_sub(t1, GR_ENTRY(ztmp, b, sz),
                            GR_ENTRY(ztmp, bstar, sz), ctx);
                    status |= _gr_dft_mul_root(GR_ENTRY(ytmp, 0, sz), t1,
                            e * kstep, inverse, rtmp, P);
                    status |= gr_add(acc, acc, GR_ENTRY(ytmp, 0, sz), ctx);
                }
            }
        }
    }

    return status;
}

/* res must not alias vec */
int
_gr_dft_mixed(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    ulong pmax = 2;
    gr_ptr ztmp, ytmp, t1, rtmp = NULL;

    if (P->n == 1)
        return gr_set(res, vec, ctx);

    for (i = 0; i < P->num_radices; i++)
        pmax = FLINT_MAX(pmax, P->radices[i]);

    ztmp = flint_malloc(2 * pmax * ctx->sizeof_elem);
    _gr_vec_init(ztmp, 2 * pmax, ctx);
    ytmp = GR_ENTRY(ztmp, pmax, ctx->sizeof_elem);

    GR_TMP_INIT(t1, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_INIT(rtmp, P->real_ctx);

    status |= _mixed(res, vec, 1, 1, P->n, 0, inverse,
            ztmp, ytmp, t1, rtmp, P, ctx);

    GR_TMP_CLEAR(t1, ctx);
    if (P->real_ctx != NULL)
        GR_TMP_CLEAR(rtmp, P->real_ctx);

    _gr_vec_clear(ztmp, 2 * pmax, ctx);
    flint_free(ztmp);

    return status;
}
