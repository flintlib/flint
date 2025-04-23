/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_leading_taylor_shift(gr_ptr shift, const gr_poly_t p,
                             const gr_poly_t q, gr_ctx_t ctx)
{
#define plc GR_ENTRY(p->coeffs, p->length - 1, ctx->sizeof_elem)
#define qlc GR_ENTRY(q->coeffs, q->length - 1, ctx->sizeof_elem)
    int status = GR_SUCCESS;

    status |= gr_zero(shift, ctx);

    if (status != GR_SUCCESS || (p == q && gr_ctx_is_canonical(ctx) == T_TRUE))
	return status;

    /* Check degree and leading coefficient */

    if (p->length > 0 && gr_is_zero(plc, ctx) != T_FALSE)
	return GR_UNABLE;
    if (q->length > 0 && gr_is_zero(qlc, ctx) != T_FALSE)
	return GR_UNABLE;

    if (p->length != q->length)
	return GR_DOMAIN;

    if (p->length == 0)
        return GR_SUCCESS;

    slong n = p->length - 1;

    status = gr_check(gr_equal(plc, qlc, ctx));
    if (status != GR_SUCCESS || n == 0)
	return status;

    /* Compute candidate shift based on coefficient of x^{q->length - 1-1} */

    gr_ptr c;

    GR_TMP_INIT(c, ctx);

    status |= gr_sub(shift, GR_ENTRY(q->coeffs, n - 1, ctx->sizeof_elem),
		            GR_ENTRY(p->coeffs, n - 1, ctx->sizeof_elem), ctx);
    status |= gr_mul_si(c, plc, n, ctx);
    status |= gr_div(shift, shift, c, ctx);

    if (status == GR_DOMAIN && gr_ctx_is_integral_domain(ctx) != T_TRUE)
        status |= GR_UNABLE;

    GR_TMP_CLEAR(c, ctx);

    return status;
#undef qlc
#undef plc
}
