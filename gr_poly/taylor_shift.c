#include "gr_poly.h"

/* todo: divconquer, convolution algorithms */
int
_gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    return _gr_poly_taylor_shift_horner(res, poly, len, c, ctx);
}

int
gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != f)
        status |= gr_poly_set(res, f, ctx);

    status |= _gr_poly_taylor_shift(res->coeffs, res->coeffs, res->length, c, ctx);
    return status;
}
