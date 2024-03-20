#include <stdio.h>
#include "gr_sparse_vec.h"

int
gr_sparse_vec_write_nz(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;
    slong nnz = vec->nnz;
    gr_stream_write(out, "[");
    for (i = 0; i < nnz; i++)
    {
        gr_stream_write(out, "\n\t");
        gr_stream_write_si(out, vec->inds[i]);
        gr_stream_write(out, ": ");
        status |= gr_write(out, GR_ENTRY(vec->nzs, i, sz), ctx);
        if (i < nnz - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "\n]");
    return status;  
}

int gr_sparse_vec_print_nz(const gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_sparse_vec_write_nz(out, vec, ctx);
}
