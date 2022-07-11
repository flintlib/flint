#include "gr_poly.h"

void
_gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz;
    truth_t eq;

    i = poly->length - 1;
    sz = ctx->sizeof_elem;

    while (i >= 0)
    {
        eq = gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);

        if (eq == T_TRUE)
        {
            GR_MUST_SUCCEED(gr_zero(GR_ENTRY(poly->coeffs, i, sz), ctx));
            i--;
        }
        else
        {
            break;
        }
    }

    poly->length = i + 1;
}
