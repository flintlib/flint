#include "gr_poly.h"

int
gr_poly_write(gr_stream_t out, const gr_poly_t poly, gr_ctx_t ctx)
{
    int status;
    slong i, n;
    slong sz;

    sz = ctx->sizeof_elem;
    n = gr_poly_length(poly, ctx);

    status = GR_SUCCESS;
    gr_stream_write(out, "[");

    for (i = 0; i < n; i++)
    {
        status |= gr_write(out, GR_ENTRY(poly->coeffs, i, sz), ctx);
        if (i < n - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return status;
}

int
gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_poly_write(out, poly, ctx);
}
