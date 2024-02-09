#include "gr_sparse_vec.h"

void
gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx)
{
    slong alloc = vec->alloc;
    /* It doesn't make sense to allocate more than the ambient dimension */
    if (nnz > vec->length)
        nnz = vec->length;
    if (nnz > alloc)
    {
        slong sz = ctx->sizeof_elem;
        if (nnz < 2 * alloc)
            nnz = 2 * alloc;
        if (nnz > vec->length)
            nnz = vec->length;
        vec->cols = flint_realloc(vec->cols, nnz * sizeof(ulong));
        vec->entries = flint_realloc(vec->entries, nnz * sz);
        _gr_vec_init(GR_ENTRY(vec->entries, alloc, sz), nnz - alloc, ctx);
        vec->alloc = len;
    }
}
