#include "gr_sparse_vec.h"

void
gr_sparse_vec_shrink_to_nnz(gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    slong nnz = vec->nnz;
    slong sz = ctx->sizeof_elem;
    if (vec->alloc > nnz)
    {
        vec->inds = flint_realloc(vec->inds, nnz * sizeof(ulong));
        _gr_vec_clear(GR_ENTRY(vec->nzs, nnz, sz), vec->alloc - nnz, ctx);
        vec->nzs = flint_realloc(vec->nzs, nnz * sz);
        vec->alloc = nnz;
    }
}
