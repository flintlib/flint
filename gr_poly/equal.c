#include "gr_poly.h"

int
_gr_poly_equal(int * res, gr_srcptr poly1, slong len1,
        gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    int status, res2, status2;
    slong sz = ctx->sizeof_elem;

    status = _gr_vec_equal(res, poly1, poly2, len2, ctx);

    if (len1 == len2 || (status == GR_SUCCESS && *res == 0))
        return status;

    status2 = _gr_vec_is_zero(&res2, GR_ENTRY(poly1, len2, sz), len2 - len1, ctx);

    if (status2 == GR_SUCCESS && res2 == 0)
    {
        *res = 0;
        return GR_SUCCESS;
    }

    if (status2 == GR_SUCCESS && res2 == 1 && status == GR_SUCCESS && *res == 0)
    {
        *res = 1;
        return GR_SUCCESS;
    }

    *res = 0;
    return GR_UNABLE;
}

int
gr_poly_equal(int * res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;

    if (len1 >= len2)
        return _gr_poly_equal(res, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    else
        return _gr_poly_equal(res, poly2->coeffs, len2, poly1->coeffs, len1, ctx);
}
