/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "gr_special.h"
#include "gr_vec.h"

/* todo: overloads */

int gr_dirichlet_chi_fmpz(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, const fmpz_t n, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_CC_ACB)
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        acb_dirichlet_chi(res, G, chi, fmpz_fdiv_ui(n, G->q), prec);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int gr_dirichlet_chi_vec(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_CC_ACB)
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_dirichlet_chi_vec(res, G, chi, len, prec);
        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int gr_dirichlet_l(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr s, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_CC_ACB)
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_dirichlet_l(res, s, G, chi, prec);
        return GR_SUCCESS;
    }

    if (ctx->which_ring == GR_CTX_RR_ARB)
    {
        acb_t t;
        slong prec;
        int status;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_init(t);
        acb_set_arb(t, s);
        acb_dirichlet_l(t, t, G, chi, prec);
        if (arb_is_zero(acb_imagref(t)) && arb_is_finite(acb_realref(t)))
        {
            arb_swap(res, acb_realref(t));
            status = GR_SUCCESS;
        }
        else
        {
            status = GR_UNABLE;
        }
        acb_clear(t);
        return status;
    }

    return GR_UNABLE;
}

int gr_dirichlet_hardy_z(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr s, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_CC_ACB)
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_dirichlet_hardy_z(res, s, G, chi, 1, prec);
        return GR_SUCCESS;
    }

    if (ctx->which_ring == GR_CTX_RR_ARB)
    {
        acb_t t;
        slong prec;
        int status;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_init(t);
        acb_set_arb(t, s);
        acb_dirichlet_hardy_z(t, t, G, chi, 1, prec);
        if (arb_is_zero(acb_imagref(t)) && arb_is_finite(acb_realref(t)))
        {
            arb_swap(res, acb_realref(t));
            status = GR_SUCCESS;
        }
        else
        {
            status = GR_UNABLE;
        }
        acb_clear(t);
        return status;
    }

    return GR_UNABLE;
}

int gr_dirichlet_hardy_theta(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr s, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_CC_ACB)
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_dirichlet_hardy_theta(res, s, G, chi, 1, prec);
        return GR_SUCCESS;
    }

    if (ctx->which_ring == GR_CTX_RR_ARB)
    {
        acb_t t;
        slong prec;
        int status;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        acb_init(t);
        acb_set_arb(t, s);
        acb_dirichlet_hardy_theta(t, t, G, chi, 1, prec);
        if (arb_is_zero(acb_imagref(t)) && arb_is_finite(acb_realref(t)))
        {
            arb_swap(res, acb_realref(t));
            status = GR_SUCCESS;
        }
        else
        {
            status = GR_UNABLE;
        }
        acb_clear(t);
        return status;
    }

    return GR_UNABLE;
}
