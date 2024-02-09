#include "gr_sparse_vec.h"

void
gr_sparse_vec_shrink_to_nnz(gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    slong nnz = vec->nnz;
    slong sz = ctx->sizeof_elem;
    if (vec->alloc > nnz)
    {
        vec->cols = flint_realloc(vec->cols, nnz * sizeof(ulong));
        _gr_vec_clear(GR_ENTRY(vec->entries, nnz, sz), vec->alloc - nnz, ctx);
        vec->entries = flint_realloc(vec->entries, nnz * sz);
        vec->alloc = nnz;
    }
}