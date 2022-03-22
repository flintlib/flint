#include "gr_vec.h"

int gr_vec_init(gr_vec_t vec, slong len, gr_ctx_t ctx)
{
    vec->length = vec->alloc = len;

    if (len == 0)
    {
        vec->entries = NULL;
    }
    else
    {
        slong sz = ctx->sizeof_elem;
        vec->entries = flint_malloc(len * sz);
        _gr_vec_init(vec->entries, len, ctx);
    }

    return GR_SUCCESS;
}
