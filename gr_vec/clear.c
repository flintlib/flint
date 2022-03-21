#include "gr_vec.h"

int gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx)
{
    if (vec->alloc != 0)
    {
        _gr_vec_clear(vec->entries, vec->alloc, ctx);
        flint_free(vec->entries);
    }

    return GR_SUCCESS;
}
