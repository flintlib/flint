#include "gr_poly.h"

int
gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    gr_poly_fit_length(poly, len, ctx);
    _gr_vec_randtest(poly->coeffs, state, len, ctx);
    _gr_poly_set_length(poly, len, ctx);
    _gr_poly_normalise(poly, ctx);
    return GR_SUCCESS;
}
