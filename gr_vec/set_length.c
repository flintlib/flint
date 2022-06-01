#include "gr_vec.h"

void
gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx)
{
    int status;

    if (vec->length > len)
    {
        status = _gr_vec_zero(GR_ENTRY(vec->entries, len, ctx->sizeof_elem), vec->length - len, ctx);
    }
    else if (vec->length < len)
    {
        gr_vec_fit_length(vec, len, ctx);
        status = _gr_vec_zero(GR_ENTRY(vec->entries, vec->length, ctx->sizeof_elem), len - vec->length, ctx);
    }

    /* todo: must handle */
    FLINT_ASSERT(status == GR_SUCCESS)

    vec->length = len;
}
