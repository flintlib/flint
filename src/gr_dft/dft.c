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

/* NOTE: this dispatcher (and the raw entry below) must route the
   Cooley-Tukey and split-radix algorithms through their _threaded
   entries; the plain _gr_dft_ct and _gr_dft_split are the serial
   kernels, also used for the sub-transforms inside worker threads,
   and never acquire threads themselves. */
static int
_gr_dft_dispatch(gr_ptr res, gr_srcptr vec, int inverse, int scrambled,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    switch (P->alg)
    {
        case GR_DFT_ALG_SPLIT:
            return _gr_dft_split_threaded(res, vec, inverse, P, ctx);
        case GR_DFT_ALG_MIXED:
            return _gr_dft_mixed(res, vec, inverse, P, ctx);
        case GR_DFT_ALG_PFA:
            return _gr_dft_pfa(res, vec, inverse, P, ctx);
        case GR_DFT_ALG_BLUESTEIN:
            return _gr_dft_bluestein(res, vec, inverse, P, ctx);
        default:
            return _gr_dft_naive(res, vec, inverse, scrambled, P, ctx);
    }
}

/* Apply the plan without the 1/n scaling of the inverse transform. */
int
_gr_dft_precomp_raw(gr_ptr res, gr_srcptr vec, int inverse,
        const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    int scrambled = (P->flags & GR_DFT_SCRAMBLED) != 0;

    if (P->alg == GR_DFT_ALG_CT || P->alg == GR_DFT_ALG_BAILEY)
    {
        /* in-place algorithms */
        if (res != vec)
            status |= _gr_vec_set(res, vec, P->n, ctx);

        if (P->alg == GR_DFT_ALG_BAILEY)
            status |= _gr_dft_bailey(res, inverse, scrambled, P, ctx);
        else
            status |= _gr_dft_ct_threaded(res, 1, inverse, scrambled, P, ctx);
    }
    else
    {
        /* out-of-place algorithms */
        if (res != vec)
        {
            status |= _gr_dft_dispatch(res, vec, inverse, scrambled, P, ctx);
        }
        else
        {
            gr_ptr tmp;
            tmp = flint_malloc(P->n * ctx->sizeof_elem);
            _gr_vec_init(tmp, P->n, ctx);
            status |= _gr_dft_dispatch(tmp, vec, inverse, scrambled, P, ctx);
            status |= _gr_vec_set(res, tmp, P->n, ctx);
            _gr_vec_clear(tmp, P->n, ctx);
            flint_free(tmp);
        }
    }

    return status;
}

int
gr_dft_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    return _gr_dft_precomp_raw(res, vec, 0, P, ctx);
}

int
gr_dft_inverse_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= _gr_dft_precomp_raw(res, vec, 1, P, ctx);

    /* scale by 1/n */
    if (P->n >= 2)
        status |= _gr_vec_div_scalar_ui(res, res, P->n, P->n, ctx);

    return status;
}

int
gr_dft(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_dft_pre_t P;

    status |= gr_dft_precomp_init(P, n, GR_DFT_ALG_AUTO, 0, ctx);

    if (status == GR_SUCCESS)
        status |= gr_dft_precomp(res, vec, P, ctx);

    gr_dft_precomp_clear(P);
    return status;
}

int
gr_dft_inverse(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_dft_pre_t P;

    status |= gr_dft_precomp_init(P, n, GR_DFT_ALG_AUTO, 0, ctx);

    if (status == GR_SUCCESS)
        status |= gr_dft_inverse_precomp(res, vec, P, ctx);

    gr_dft_precomp_clear(P);
    return status;
}
