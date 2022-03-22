#include "gr_poly.h"

int
_gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz;
    int status;
    truth_t eq;

    i = poly->length - 1;
    sz = ctx->sizeof_elem;

    status = GR_SUCCESS;

    while (i >= 0)
    {
        eq = gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);

        if (eq == T_TRUE)
        {
            gr_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);
            i--;
        }
        else
        {
            if (eq == T_UNKNOWN)
                status = GR_UNABLE;

            break;
        }
    }

    poly->length = i + 1;
    return status;
}
