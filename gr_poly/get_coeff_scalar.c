#include "gr_poly.h"

int
gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong i, gr_ctx_t ctx)
{
    if (i < 0 || i >= poly->length)
        return gr_zero(res, ctx);
    else
        return gr_set(res, GR_ENTRY(poly->coeffs, i, ctx->sizeof_elem), ctx);
}
