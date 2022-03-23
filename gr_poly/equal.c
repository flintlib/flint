#include "gr_poly.h"

truth_t
_gr_poly_equal(gr_srcptr poly1, slong len1,
        gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    truth_t eq, eq2;
    slong sz = ctx->sizeof_elem;

    eq = _gr_vec_equal(poly1, poly2, len2, ctx);
    if (len1 == len2 || eq == T_FALSE)
        return eq;

    eq2 = _gr_vec_is_zero(GR_ENTRY(poly1, len2, sz), len1 - len2, ctx);

    return truth_and(eq, eq2);
}

truth_t
gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;

    if (len1 >= len2)
        return _gr_poly_equal(poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    else
        return _gr_poly_equal(poly2->coeffs, len2, poly1->coeffs, len1, ctx);
}
