#include "gr_poly.h"

int
_gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz;
    int status, is_zero;

    i = poly->length - 1;
    sz = ctx->sizeof_elem;

    status = GR_SUCCESS;

    while (i >= 0)
    {
        status = gr_is_zero(&is_zero, GR_ENTRY(poly->coeffs, i, sz), ctx);

        if (status == GR_SUCCESS && is_zero)
        {
            gr_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);
            i--;
        }
        else
        {
            break;
        }
    }

    poly->length = i + 1;
    return status;
}
