#include <stdlib.h>
#include "gr_sparse_vec.h"


static int gr_sparse_vec_ulong_cmp(const void* a, const void* b)
{
    ulong av = *((ulong*)(a));
    ulong bv = *((ulong*)(b));
    return (av < bv ? -1 : (av > bv ? 1 : 0));
}


gr_ptr
gr_sparse_vec_find_entry(gr_sparse_vec_t vec, slong col, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    ulong* bs = NULL;
    if (col < 0 || col >= vec->length)
        return NULL;
    bs = bsearch(&col, vec->inds, vec->nnz, sizeof(slong), gr_sparse_vec_ulong_cmp);
    if (bs == NULL)
        return NULL;
    return GR_ENTRY(vec->nzs, bs - vec->inds, sz);
}

int
gr_sparse_vec_get_entry(gr_ptr res, gr_sparse_vec_t vec, slong col, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    ulong* bs = NULL;
    if (col < 0 || col >= vec->length)
        return GR_DOMAIN;
    bs = bsearch(&col, vec->inds, vec->nnz, sizeof(slong), gr_sparse_vec_ulong_cmp);
    if (bs == NULL)
        return gr_zero(res, ctx);
    return gr_set(res, GR_ENTRY(vec->nzs, bs - vec->inds, sz), ctx);
}
