#include "gr_poly.h"

int
gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
{
    if (res != src)
    {
        gr_poly_fit_length(res, src->length, ctx);
        _gr_vec_set(res->coeffs, src->coeffs, src->length, ctx);
        _gr_poly_set_length(res, src->length, ctx);
    }

    return GR_SUCCESS;
}
