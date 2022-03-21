#include "gr_poly.h"

int
gr_poly_neg_one(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_fit_length(res, 1, ctx);
    _gr_poly_set_length(res, 1, ctx);
    return gr_neg_one(res->coeffs, ctx);
}
