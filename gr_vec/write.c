#include "gr_vec.h"

int
_gr_vec_write(gr_stream_t out, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;

    gr_stream_write(out, "[");
    for (i = 0; i < len; i++)
    {
        status |= gr_write(out, GR_ENTRY(vec, i, sz), ctx);
        if (i < len - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return GR_SUCCESS;
}

int
gr_vec_write(gr_stream_t out, const gr_vec_t vec, gr_ctx_t ctx)
{
    return _gr_vec_write(out, vec->entries, vec->length, ctx);
}

int
gr_vec_print(const gr_vec_t vec, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_vec_write(out, vec, ctx);
}
