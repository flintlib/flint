#include "gr_sparse_vec.h"

void
gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx)
{
    slong alloc = vec->alloc;
    slong new_alloc = nnz;
    /* It doesn't make sense to allocate more than the ambient dimension */
    if (new_alloc > vec->length)
        new_alloc = vec->length;
    if (new_alloc > alloc)
    {
        slong sz = ctx->sizeof_elem;
        if (new_alloc < 2 * alloc)
            new_alloc = 2 * alloc;
        if (new_alloc > vec->length)
            new_alloc = vec->length;
        vec->cols = flint_realloc(vec->cols, new_alloc * sizeof(ulong));
        vec->entries = flint_realloc(vec->entries, new_alloc * sz);
        _gr_vec_init(GR_ENTRY(vec->entries, alloc, sz), new_alloc - alloc, ctx);
        vec->alloc = new_alloc;
    }
}
