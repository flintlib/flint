#include "gr_poly.h"

int
gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)
{
    slong alloc = poly->alloc;

    if (len > alloc)
    {
        slong sz = ctx->sizeof_elem;

        if (len < 2 * alloc)
            len = 2 * alloc;

        poly->coeffs = flint_realloc(poly->coeffs, len * sz);
        _gr_vec_init(GR_ENTRY(poly->coeffs, poly->alloc, sz), len - alloc, ctx);
        poly->alloc = len;
    }

    return GR_SUCCESS;
}
