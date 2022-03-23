#include "gr_vec.h"

void gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx)
{
    if (vec->alloc != 0)
    {
        _gr_vec_clear(vec->entries, vec->alloc, ctx);
        flint_free(vec->entries);
    }
}
