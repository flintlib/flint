/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_poly.h"
#include "gr_vec.h"

int
gr_generic_poly_factor_roots(
        gr_ptr c, gr_vec_t fac, gr_vec_t mult, gr_srcptr elt, int flags,
        gr_ctx_t ctx)
{
    gr_ctx_t pctx;
    gr_vec_t roots;
    fmpz_t tot;

    gr_poly_struct * poly = (gr_poly_struct *) elt;

    if (gr_ctx_is_integral_domain(ctx) != T_TRUE)
        return GR_UNABLE;

    if (poly->length == 0)
        return GR_DOMAIN;

    gr_ctx_init_gr_poly(pctx, ctx);
    gr_vec_init(roots, 0, ctx);
    fmpz_init(tot);

    int status = GR_SUCCESS;

    gr_ptr lc = gr_poly_coeff_ptr(poly, poly->length - 1, ctx);
    if (gr_is_zero(lc, ctx) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }
    status |= gr_poly_set_scalar(c, lc, ctx);

    status |= gr_poly_roots(roots, mult, poly, 0, ctx);

    gr_ctx_t ZZ;
    gr_ctx_init_fmpz(ZZ);
    GR_MUST_SUCCEED(_gr_vec_sum(tot, mult->entries, mult->length, ZZ));
    gr_ctx_clear(ZZ);

    if (!fmpz_equal_si(tot, poly->length - 1))
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    gr_vec_set_length(fac, roots->length, pctx);
    for (slong i = 0; i < roots->length; i++)
    {
        gr_ptr rt = gr_vec_entry_ptr(roots, i, ctx);
        gr_poly_struct * f = gr_vec_entry_ptr(fac, i, pctx);
        status |= gr_poly_gen(f, ctx);
        status |= gr_neg(rt, rt, ctx);
        gr_swap(f->coeffs, rt, ctx);
    }

cleanup:

    fmpz_clear(tot);
    gr_vec_clear(roots, ctx);
    gr_ctx_clear(pctx);

    return status;
}

