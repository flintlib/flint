#include "gr_poly.h"

void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx)
{
    gr_poly_init(poly, ctx);
    gr_poly_fit_length(poly, len, ctx);
}
