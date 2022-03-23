#include "gr_vec.h"

void
gr_vec_fit_length(gr_vec_t vec, slong len, gr_ctx_t ctx)
{
    slong alloc = vec->alloc;

    if (len > alloc)
    {
        slong sz = ctx->sizeof_elem;

        if (len < 2 * alloc)
            len = 2 * alloc;

        vec->entries = flint_realloc(vec->entries, len * sz);
        _gr_vec_init(GR_ENTRY(vec->entries, vec->alloc, sz), len - alloc, ctx);
        vec->alloc = len;
    }
}
