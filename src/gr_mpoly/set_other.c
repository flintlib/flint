/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpoly.h"
#include "fmpz_mpoly.h"
#include "gr_mpoly.h"

static int
_same_vars(char ** v1, char ** v2, slong n)
{
    slong i;

    if (v1 == NULL && v2 == NULL)
        return 1;

    if (v1 == NULL || v2 == NULL)
        return 0;

    for (i = 0; i < n; i++)
        if (strcmp(v1[i], v2[i]))
            return 0;

    return 1;
}

int
_gr_mpoly_set_gr_mpoly_other(gr_mpoly_t res, const gr_mpoly_t A, gr_mpoly_ctx_t A_ctx, gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);

    mpoly_ctx_struct * A_mctx = GR_MPOLY_MCTX(A_ctx);
    gr_ctx_struct * A_cctx = GR_MPOLY_CCTX(A_ctx);

    slong A_nvars, nvars;
    slong len;

    len = A->length;
    nvars = GR_MPOLY_NVARS(ctx);
    A_nvars = GR_MPOLY_NVARS(A_ctx);

    if (A->length == 0)
    {
        /* Before converting the zero polynomial to the zero polynomial,
           make sure that 0 -> 0 is a legal conversion between the
           base rings. */
        if (A_cctx == cctx)
        {
            return gr_mpoly_zero(res, ctx);
        }
        else
        {
            int status = GR_SUCCESS;
            gr_ptr c, d;

            GR_TMP_INIT(c, A_cctx);
            GR_TMP_INIT(d, cctx);

            status |= gr_mpoly_zero(res, ctx);
            status |= gr_set_other(d, c, A_cctx, cctx);

            GR_TMP_CLEAR(c, A_cctx);
            GR_TMP_CLEAR(d, cctx);
            return status;
        }
    }

    /* Simplest case: same variables & term ordering */
    if (nvars == A_nvars
        && mctx->ord == A_mctx->ord
        && _same_vars(GR_MPOLY_VARS(A_ctx), GR_MPOLY_VARS(ctx), nvars))
    {
        int status = GR_SUCCESS;
        slong i, N, sz = cctx->sizeof_elem, A_sz = A_cctx->sizeof_elem;

        N = mpoly_words_per_exp(A->bits, mctx);
        gr_mpoly_fit_length_reset_bits(res, len, A->bits, ctx);

        for (i = 0; status == GR_SUCCESS && i < len; i++)
        {
            status |= gr_set_other(GR_ENTRY(res->coeffs, i, sz), GR_ENTRY(A->coeffs, i, A_sz), A_cctx, cctx);
        }

        if (status == GR_SUCCESS)
        {
            mpoly_copy_monomials(res->exps, A->exps, len, N);
            _gr_mpoly_set_length(res, len, ctx);

            /* there may be zero coefficients; remove them
               (todo: do inline in first loop) */
            /* if (cctx != A_cctx) */
            status = gr_mpoly_combine_like_terms(res, ctx);
        }
        else
        {
            _gr_mpoly_set_length(res, 0, ctx);
        }

        return status;
    }

    /* Handle the case where the used variables in A are a subset of the variables in ctx */
    {
        int * used_mask;
        slong len, num_used;
        slong * used;
        slong * translation;
        int status = GR_SUCCESS;
        slong i, j, N, sz = cctx->sizeof_elem, A_sz = A_cctx->sizeof_elem;
        slong num_translated;

        len = A->length;

        /* Both contexts have named variables */
        if (GR_MPOLY_VARS(ctx) != NULL && GR_MPOLY_VARS(A_ctx) != NULL)
        {
            used_mask = flint_calloc(A_nvars, sizeof(int));
            used = flint_malloc(sizeof(slong) * A_nvars);
            translation = flint_malloc(sizeof(slong) * A_nvars);

            mpoly_used_vars_or(used_mask, A->exps, len, A->bits, A_mctx);

            num_used = 0;
            for (i = 0; i < A_nvars; i++)
            {
                if (used_mask[i] != 0)
                {
                    used[num_used] = i;
                    num_used++;
                }
            }

            if (num_used <= nvars)
            {
                /* TODO: use a non-quadratic match algorithm */
                num_translated = 0;

                for (i = 0; i < num_used; i++)
                {
                    for (j = 0; j < nvars; j++)
                    {
                        if (!strcmp(GR_MPOLY_VARS(A_ctx)[used[i]], GR_MPOLY_VARS(ctx)[j]))
                        {
                            translation[num_translated] = j;
                            num_translated++;
                            break;
                        }
                    }
                }

                if (num_translated != num_used)
                {
                    /* todo: GR_DOMAIN in appropriate cases */
                    status = GR_UNABLE;
                }
                else if (A->bits <= FLINT_BITS)
                {
                    ulong * A_exps;
                    ulong * exps;

                    A_exps = flint_malloc(sizeof(ulong) * A_nvars);
                    exps = flint_calloc(nvars, sizeof(ulong));

                    N = mpoly_words_per_exp(A->bits, A_mctx);
                    gr_mpoly_fit_length_reset_bits(res, len, A->bits, ctx);

                    for (i = 0; status == GR_SUCCESS && i < len; i++)
                    {
                        status |= gr_set_other(GR_ENTRY(res->coeffs, i, sz), GR_ENTRY(A->coeffs, i, A_sz), A_cctx, cctx);
                        mpoly_get_monomial_ui(A_exps, A->exps + i * N, A->bits, A_mctx);
                        for (j = 0; j < num_used; j++)
                            exps[translation[j]] = A_exps[used[j]];
                        _gr_mpoly_push_exp_ui(res, exps, ctx);
                    }

                    _gr_mpoly_set_length(res, len, ctx);

                    /* term order may be different */
                    gr_mpoly_sort_terms(res, ctx);
                    /* todo: combine zero checks with main loop */
                    status |= gr_mpoly_combine_like_terms(res, ctx);

                    flint_free(A_exps);
                    flint_free(exps);
                }
                else
                {
                    /* todo: fmpz exponents */
                    status = GR_UNABLE;
                }
            }

            flint_free(translation);
            flint_free(used_mask);
            flint_free(used);
        }
        else
        {
            status = GR_UNABLE;
        }

        return status;
    }

    /* Other cases are not yet handled, e.g. R[x,y][s,t] -> R[x,y,s,t] */

    return GR_UNABLE;
}


/* hack: to convert an fmpz_mpoly, mock up a gr_mpoly over ZZ */

typedef struct
{
    fmpz_mpoly_ctx_t mctx;
    char ** vars;
}
_gr_fmpz_mpoly_ctx_t;

#define _FMPZ_MPOLY_CTX(ring_ctx) ((_gr_fmpz_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define _FMPZ_MPOLY_MCTX(ring_ctx) (_FMPZ_MPOLY_CTX(ring_ctx)->mctx)

int
_gr_mpoly_set_fmpz_mpoly(gr_mpoly_t res, const fmpz_mpoly_t A, gr_ctx_t A_ctx, gr_mpoly_ctx_t ctx)
{
    int status;

    gr_ctx_t ZZ;
    gr_mpoly_ctx_t A_ctx1;
    gr_mpoly_t t;

    gr_ctx_init_fmpz(ZZ);  /* no need to free */

    *A_ctx1 = *ctx;
    GR_MPOLY_MCTX(A_ctx1) = _FMPZ_MPOLY_MCTX(A_ctx)->minfo;
    GR_MPOLY_CCTX(A_ctx1) = ZZ;
    GR_MPOLY_VARS(A_ctx1) = _FMPZ_MPOLY_CTX(A_ctx)->vars;

    t->coeffs = A->coeffs;
    t->exps = A->exps;
    t->length = A->length;
    t->bits = A->bits;
    t->coeffs_alloc = 0;
    t->exps_alloc = 0;

    status = _gr_mpoly_set_gr_mpoly_other(res, t, A_ctx1, ctx);

    return status;
}

int gr_mpoly_set_other(gr_mpoly_t res, gr_srcptr A, gr_ctx_t A_ctx, gr_mpoly_ctx_t ctx)
{
    /* mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx); */
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);

    if (A_ctx == ctx)
    {
        return gr_mpoly_set(res, A, ctx);
    }
    else if (A_ctx == cctx)
    {
        return gr_mpoly_set_scalar(res, A, ctx);
    }
    else if (A_ctx->which_ring == GR_CTX_GR_MPOLY)
    {
        return _gr_mpoly_set_gr_mpoly_other(res, A, A_ctx, ctx);
    }
    else if (A_ctx->which_ring == GR_CTX_FMPZ_MPOLY)
    {
        return _gr_mpoly_set_fmpz_mpoly(res, A, A_ctx, ctx);
    }
    else
    {
        /* Otherwise, try to convert to the coefficient ring. */
        int status;
        gr_ptr t;

        GR_TMP_INIT(t, cctx);

        status = gr_set_other(t, A, A_ctx, cctx);

        if (status == GR_SUCCESS)
        {
            status = gr_mpoly_set_scalar(res, t, ctx);
        }
        else
        {
            /* We can't reliably return GR_DOMAIN here. */
            status = GR_UNABLE;
        }

        GR_TMP_CLEAR(t, cctx);

        return status;
    }
}
