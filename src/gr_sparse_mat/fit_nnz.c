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
