#include "gr_poly.h"

int
_gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;

    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
            gr_swap(GR_ENTRY(res, i, sz), GR_ENTRY(res, n - 1 - i, sz), ctx);

        for (i = 0; i < n - len; i++)
            gr_zero(GR_ENTRY(res, i, sz), ctx);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            gr_zero(GR_ENTRY(res, i, sz), ctx);

        for (i = 0; i < len; i++)
            gr_set(GR_ENTRY(res, (n - len) + i, sz), GR_ENTRY(poly, (len - 1) - i, sz), ctx);
    }

    return GR_SUCCESS;
}

int
gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)
{
    slong len = FLINT_MIN(n, poly->length);

    if (len == 0)
        return gr_poly_zero(res, ctx);

    gr_poly_fit_length(res, n, ctx);
    _gr_poly_reverse(res->coeffs, poly->coeffs, len, n, ctx);
    _gr_poly_set_length(res, n, ctx);
    return _gr_poly_normalise(res, ctx);
}
