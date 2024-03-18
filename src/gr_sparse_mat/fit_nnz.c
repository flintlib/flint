#include "gr_sparse_mat.h"

void
gr_csr_mat_fit_nnz(gr_csr_mat_t mat, slong nnz, gr_ctx_t ctx)
{
    slong alloc = mat->alloc;
    slong new_alloc = nnz;
    if (new_alloc > alloc)
    {
        slong sz = ctx->sizeof_elem;
        mat->cols = flint_realloc(mat->cols, new_alloc * sizeof(ulong));
        mat->entries = flint_realloc(mat->entries, new_alloc * sz);
        _gr_vec_init(GR_ENTRY(mat->entries, alloc, sz), new_alloc - alloc, ctx);
        mat->alloc = new_alloc;
    }
}

void
gr_csr_mat_shrink_to_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)
{
    slong nnz = mat->nnz;
    slong sz = ctx->sizeof_elem;
    if (mat->alloc > nnz)
    {
        mat->cols = flint_realloc(mat->cols, nnz * sizeof(ulong));
        _gr_vec_clear(GR_ENTRY(mat->entries, nnz, sz), mat->alloc - nnz, ctx);
        mat->entries = flint_realloc(mat->entries, nnz * sz);
        mat->alloc = nnz;
    }    
}

void
gr_lil_mat_shrink_to_nnz(gr_lil_mat_t mat, gr_ctx_t ctx)
{
    slong row;

    for (row = 0; row < mat->r; ++row)
    {
        gr_sparse_vec_shrink_to_nnz(mat->rows[row], ctx);
    }
}

