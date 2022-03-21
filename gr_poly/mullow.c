#include "gr_poly.h"


int
_gr_poly_mullow(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    int status;
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
        return gr_mul(res, poly1, poly2, ctx);

    if (len1 == 1)
        return _gr_vec_scalar_mul(res, poly2, n, poly1, ctx);

    if (len2 == 1)
        return _gr_vec_scalar_mul(res, poly1, n, poly2, ctx);

    /* todo: Squaring
    if (poly1 == poly2 && len1 == len2)
    {
        _gr_poly_sqrlow_classical(res, poly1, len1, n, ctx);
        return;
    } */

    /* General case */
    {
        slong i, top1, top2, sz;
        sz = ctx->sizeof_elem;

        status = gr_mul(res, poly1, poly2, ctx);

        for (i = 1; i < n; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), NULL, 0, GR_ENTRY(poly1, i - top2, sz), GR_ENTRY(poly2, i - top1, sz), top1 + top2 - i + 1, ctx);
        }

        return status;
    }
}

int
gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong n, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;
    n = FLINT_MIN(n, len_out);

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status = _gr_poly_mullow(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
