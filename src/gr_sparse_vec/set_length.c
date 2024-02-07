#include "gr_sparse_vec.h"

void
gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx)
{
    slong i;
    vec->len = len;
    for (i=vec->nnz-1; i>=0; i--)
    {
        
    }

    slong alloc = vec->alloc;
    if (len > alloc)
    {
        slong sz = ctx->sizeof_elem;
        if (len < 2 * alloc)
            len = 2 * alloc;
        vec->inds = flint_realloc(vec->inds, len * sizeof(ulong));
        vec->entries = flint_realloc(vec->entries, len * sz);
        _gr_vec_init(GR_ENTRY(vec->entries, vec->alloc, sz), len - alloc, ctx);
        vec->alloc = len;
    }
}



#include "gr_vec.h"

void
gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx)
{
    if (vec->length > len)
    {
        /* potentially free no longer used elements */
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(vec->entries, len, ctx->sizeof_elem), vec->length - len, ctx));
    }
    else if (vec->length < len)
    {
        gr_vec_fit_length(vec, len, ctx);
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(vec->entries, vec->length, ctx->sizeof_elem), len - vec->length, ctx));
    }

    vec->length = len;
}
