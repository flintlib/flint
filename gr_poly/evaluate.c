#include "gr_poly.h"

int
_gr_poly_evaluate(gr_ptr y, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_horner(y, f, len, x, ctx);
}

int
gr_poly_evaluate(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t ctx)
{
    return _gr_poly_evaluate(res, f->coeffs, f->length, a, ctx);
}
