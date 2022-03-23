#include "gr_poly.h"

void
_gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx)
{
    if (poly->length > len)
    {
        _gr_vec_zero(GR_ENTRY(poly->coeffs, len, ctx->sizeof_elem), poly->length - len, ctx);
    }

    poly->length = len;
}
