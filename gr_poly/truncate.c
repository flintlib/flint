#include "gr_poly.h"

int
gr_poly_truncate(gr_poly_t poly, slong newlen, gr_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        _gr_poly_set_length(poly, newlen, ctx);
        _gr_poly_normalise(poly, ctx);
    }

    return GR_SUCCESS;
}
