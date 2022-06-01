#include "gr_poly.h"

void
_gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz;
    truth_t eq;
    int status;

    i = poly->length - 1;
    sz = ctx->sizeof_elem;

    while (i >= 0)
    {
        eq = gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);

        if (eq == T_TRUE)
        {
            status = gr_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);
            /* todo: must handle status (if zero can fail) */
            status = status;
            i--;
        }
        else
        {
            break;
        }
    }

    poly->length = i + 1;
}
