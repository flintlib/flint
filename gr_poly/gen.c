#include "gr_poly.h"

int
gr_poly_gen(gr_poly_t res, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_poly_zero(res, ctx);
    status |= gr_poly_set_coeff_ui(res, 1, 1, ctx);
    return status;
}
