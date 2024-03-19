#include <stdio.h>
#include "gr_sparse_mat.h"

int
gr_csr_mat_write_nz(gr_stream_t out, const gr_csr_mat_t mat, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong row;
    gr_sparse_vec_t tmp;
    gr_stream_write(out, "[");
    for (row = 0; row < mat->r; row++)
    {
        _gr_csr_mat_borrow_row(tmp, mat, row, ctx);
        status |= gr_sparse_vec_write_nz(out, tmp, ctx);
        if (row < mat->r - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return status;  
}

int
gr_lil_mat_write_nz(gr_stream_t out, const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong row;

    gr_stream_write(out, "[");
    for (row = 0; row < mat->r; row++)
    {
        status |= gr_sparse_vec_write_nz(out, &mat->rows[row], ctx);
        if (row < mat->r - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return status;  
}

int
gr_coo_mat_write_nz(gr_stream_t out, const gr_coo_mat_t mat, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong nz;
    slong sz = ctx->sizeof_elem;

    gr_stream_write(out, "[");
    for (nz = 0; nz < mat->nnz; nz++)
    {
        gr_stream_write(out, "(");
        gr_stream_write_si(out, mat->rows[nz]);
        gr_stream_write(out, ", ");
        gr_stream_write_si(out, mat->cols[nz]);
        gr_stream_write(out, "): ");
        status |= gr_write(out, GR_ENTRY(mat->nzs, nz, sz), ctx);
        if (nz < mat->nnz - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return status;  
}

int gr_csr_mat_print_nz(const gr_csr_mat_t mat, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_csr_mat_write_nz(out, mat, ctx);
}

int gr_lil_mat_print_nz(const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_lil_mat_write_nz(out, mat, ctx);
}

int gr_coo_mat_print_nz(const gr_coo_mat_t mat, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_coo_mat_write_nz(out, mat, ctx);
}


