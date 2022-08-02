#include "gr_poly.h"

int
_gr_poly_evaluate_other(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_other_horner(res, f, len, x, x_ctx, ctx);
}

int
gr_poly_evaluate_other(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t a_ctx, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_other(res, f->coeffs, f->length, a, a_ctx, ctx);
}
