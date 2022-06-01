#include "gr_vec.h"

int
gr_vec_set(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != src)
    {
        gr_vec_set_length(res, src->length, ctx);
        status = _gr_vec_set(res->entries, src->entries, src->length, ctx);
    }

    return status;
}
