/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

int gr_mpoly_mul(gr_mpoly_t poly1,
    const gr_mpoly_t poly2,
    const gr_mpoly_t poly3,
    gr_mpoly_ctx_t ctx)
{
    slong len2 = poly2->length;
    slong len3 = poly3->length;

    if (len2 == 0 || len3 == 0)
        return gr_mpoly_zero(poly1, ctx);

    if (len3 == 1)
        return gr_mpoly_mul_monomial(poly1, poly2, poly3, ctx);

    /* todo: could have a version of mul_monomial for the noncommutative case */
    if (len2 == 1 && gr_ctx_is_approx_commutative_ring(GR_MPOLY_CCTX(ctx)) == T_TRUE)
        return gr_mpoly_mul_monomial(poly1, poly3, poly2, ctx);

    if (len2 * len3 > ctx->size_limit)
        return GR_UNABLE | gr_mpoly_zero(poly1, ctx);

    if (len2 * len3 > 10000 && len2 > 2 && len3 > 2)
        return gr_mpoly_mul_heap_threaded(poly1, poly2, poly3, ctx);
    else
        return gr_mpoly_mul_johnson(poly1, poly2, poly3, ctx);
}
