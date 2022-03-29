#include "gr_poly.h"

/* todo: move out */
/* todo: alternative algorithm */
int
_gr_vec_set_powers(gr_ptr xs, gr_srcptr x, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz = ctx->sizeof_elem;

    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            status |= gr_one(GR_ENTRY(xs, i, sz), ctx);
        else if (i == 1)
            status |= gr_set(GR_ENTRY(xs, i, sz), x, ctx);
        else if (i % 2 == 0)
            status |= gr_sqr(GR_ENTRY(xs, i, sz), GR_ENTRY(xs, i / 2, sz), ctx);
        else
            status |= gr_mul(GR_ENTRY(xs, i, sz), GR_ENTRY(xs, i - 1, sz), x, ctx);
    }

    return status;
}

int
_gr_poly_evaluate_rectangular(gr_ptr y, gr_srcptr poly,
    slong len, gr_srcptr x, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len <= 2)
    {
        if (len == 0)
            return gr_zero(y, ctx);

        if (len == 1)
            return gr_set(y, poly, ctx);

        status |= gr_mul(y, x, GR_ENTRY(poly, 1, sz), ctx);
        status |= gr_add(y, y, GR_ENTRY(poly, 0, sz), ctx);
        return status;
    }
    else
    {
        slong i, m, r;
        gr_ptr xs;
        gr_ptr s, t, c;
        GR_TMP_START;

        m = n_sqrt(len) + 1;
        r = (len + m - 1) / m;

        GR_TMP_INIT_VEC(xs, m + 1, ctx);
        GR_TMP_INIT3(s, t, c, ctx);

        status |= _gr_vec_set_powers(xs, x, m + 1, ctx);
        status |= _gr_vec_dot(y, GR_ENTRY(poly, (r - 1) * m, sz), 0, GR_ENTRY(xs, 1, sz),
            GR_ENTRY(poly, (r - 1) * m + 1, sz), len - (r - 1) * m - 1, ctx);

        for (i = r - 2; i >= 0; i--)
        {
            status |= _gr_vec_dot(s, GR_ENTRY(poly, i * m, sz), 0, GR_ENTRY(xs, 1, sz),
                GR_ENTRY(poly, i * m + 1, sz), m - 1, ctx);
            status |= gr_mul(y, y, GR_ENTRY(xs, m, sz), ctx);
            status |= gr_add(y, y, s, ctx);
        }

        GR_TMP_CLEAR_VEC(xs, m + 1, ctx);
        GR_TMP_CLEAR3(s, t, c, ctx);

        GR_TMP_END;
        return status;
    }
}

int
gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_rectangular(res, f->coeffs, f->length, x, ctx);
}
