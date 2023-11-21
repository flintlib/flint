/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"
#include "partitions.h"
#include "arb.h"
#include "gr_special.h"

int
gr_generic_partitions_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        partitions_fmpz_ui(res, n);
        return GR_SUCCESS;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_partitions_ui(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else    /* compute via fmpz */
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        partitions_fmpz_ui(t, n);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_partitions_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        partitions_fmpz_fmpz(res, n, 0);
        return GR_SUCCESS;
    }
    else if (fmpz_sgn(n) < 0)
    {
        return gr_zero(res, ctx);
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_partitions_fmpz(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else    /* compute via fmpz */
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        partitions_fmpz_fmpz(t, n, 0);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

/* xxx */
#define NMOD_CTX_REF(ring_ctx) (((nmod_t *)((ring_ctx))))
#define NMOD_CTX(ring_ctx) (*NMOD_CTX_REF(ring_ctx))

int
gr_generic_partitions_vec(gr_ptr res, slong len, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_number_of_partitions_vec(res, len);
        return GR_SUCCESS;
    }

    if (ctx->which_ring == GR_CTX_NMOD)
    {
        arith_number_of_partitions_nmod_vec(res, len, NMOD_CTX(ctx));
        return GR_SUCCESS;
    }

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute numerically via arb */
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        if (len > prec)     /* todo: better criterion */
        {
            slong i, sz = ctx->sizeof_elem;
            int status = GR_SUCCESS;
            gr_ctx_t RR;
            arb_t t;

            gr_ctx_init_real_arb(RR, prec);

            arb_init(t);
            for (i = 0; i < len; i++)
            {
                arb_partitions_ui(t, i, prec);
                status |= gr_set_other(GR_ENTRY(res, i, sz), t, RR, ctx);
            }

            arb_clear(t);
            gr_ctx_clear(RR);

            return status;
        }
    }

    /* compute exactly via Z */
    {
        int status = GR_SUCCESS;
        slong i, sz = ctx->sizeof_elem;
        gr_ctx_t ZZ;
        fmpz * t;
        gr_ctx_init_fmpz(ZZ);
        GR_TMP_INIT_VEC(t, len, ZZ);

        arith_number_of_partitions_vec(t, len);
        for (i = 0; i < len && status == GR_SUCCESS; i++)
            status |= gr_set_fmpz(GR_ENTRY(res, i, sz), t + i, ctx);

        GR_TMP_CLEAR_VEC(t, len, ZZ);
        gr_ctx_clear(ZZ);
        return status;
    }
}
