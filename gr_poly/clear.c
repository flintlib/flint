#include "gr_poly.h"

void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx)
{
    if (poly->alloc != 0)
    {
        _gr_vec_clear(poly->coeffs, poly->alloc, ctx);
        flint_free(poly->coeffs);
    }
}
