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
