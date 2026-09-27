/*
    Copyright (C) 2026 Mael Hostettler
    Copyright (C) 2026 Antoine Bak
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "n_poly.h"
#include "nmod_mpoly.h"
#include "gr_poly.h"

/* The algorithm itself is n_bpoly_mod_resultant, over the dense bivariate
   type of the nmod_mpoly module. A polynomial over a polynomial ring over
   nmod has exactly that layout: a gr_poly_struct is a pointer, an alloc and a
   length, as n_bpoly_struct is, and its coefficients are gr_poly_structs
   whose own coefficients are a flat array of ulong, as n_poly_struct is. The
   inputs are therefore passed on through a cast rather than converted. */

int
_gr_poly_resultant_multipoint_cutoff(slong lenA, slong lenB, slong npoints)
{
    return n_bpoly_mod_resultant_cutoff(lenA, lenB, npoints);
}

int
_gr_poly_resultant_multipoint(gr_ptr res, gr_srcptr A, slong lenA,
                              gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    gr_ctx_struct * cctx;

    if (ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    cctx = POLYNOMIAL_ELEM_CTX(ctx);

    if (cctx->which_ring != GR_CTX_NMOD)
        return GR_UNABLE;

    if (lenB <= 1)
        return _gr_poly_resultant_small(res, A, lenA, B, lenB, ctx);

    /* the precomputations assume a prime modulus */
    if (gr_ctx_is_field(cctx) != T_TRUE)
        return GR_UNABLE;

    FLINT_ASSERT(sizeof(gr_poly_struct) == sizeof(n_poly_struct));
    FLINT_ASSERT(cctx->sizeof_elem == sizeof(ulong));

    if (!_n_bpoly_mod_resultant((n_poly_struct *) res,
            (const n_poly_struct *) A, lenA,
            (const n_poly_struct *) B, lenB, NMOD_CTX(cctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
gr_poly_resultant_multipoint(gr_ptr r, const gr_poly_t f,
                             const gr_poly_t g, gr_ctx_t ctx)
{
    slong len1 = f->length;
    slong len2 = g->length;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len1 == 0 || len2 == 0)
    {
        return gr_zero(r, ctx);
    }

    if (gr_is_zero(GR_ENTRY(f->coeffs, len1 - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(g->coeffs, len2 - 1, sz), ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }

    if (len1 >= len2)
    {
        status |= _gr_poly_resultant_multipoint(r, f->coeffs, len1, g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_multipoint(r, g->coeffs, len2, f->coeffs, len1, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
