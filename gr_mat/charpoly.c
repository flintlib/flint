#include "gr_mat.h"

/* todo: algorithm selection */
int
_gr_mat_charpoly(gr_ptr cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    return _gr_mat_charpoly_berkowitz(cp, mat, ctx);
}

int gr_mat_charpoly(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_poly_fit_length(cp, mat->r + 1, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly(cp->coeffs, mat, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
