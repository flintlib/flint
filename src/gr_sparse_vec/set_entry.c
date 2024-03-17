#include <stdlib.h>
#include "gr_sparse_vec.h"

static int gr_sparse_vec_ulong_cmp(const void* a, const void* b)
{
    ulong av = *((ulong*)(a));
    ulong bv = *((ulong*)(b));
    return (av < bv ? -1 : (av > bv ? 1 : 0));
}

int
gr_sparse_vec_set_entry(gr_sparse_vec_t vec, slong col, gr_srcptr entry, gr_ctx_t ctx)
{
    slong i,j;
    slong sz = ctx->sizeof_elem;
    slong nnz = vec->nnz;
    ulong* bs = NULL;
    if (col < 0 || col >= vec->length)
        return GR_DOMAIN;
    bs = bsearch(&col, vec->inds, vec->nnz, sizeof(slong), gr_sparse_vec_ulong_cmp);
    if (bs != NULL)
    {
        i = (bs - vec->inds);
        return gr_set(GR_ENTRY(vec->entries, i, sz), entry, ctx);
    }
    else
    {
        /* In this case, it's going to cost linear time anyway, so just brute force search */
        gr_sparse_vec_fit_nnz(vec, vec->nnz+1, ctx);
        for (i = 0; i < nnz; i++)
        {
            if (vec->inds[i] < col)
                break;
        }
        /* Now if i < nnz, then we should put the entry right *after* it */
        /* If i == nnz, we should put the entry right there */
        if (i < nnz)
            i++;
        if (i < nnz)
        {
            /* If necessary, move everything up to make space */
            memmove(vec->inds + i + 1, vec->inds + i, (vec->nnz - i)*sizeof(slong));
            for (j = vec->nnz; j > i; j--)
            {
                gr_swap(GR_ENTRY(vec->entries, j-1, sz), GR_ENTRY(vec->entries, j, sz), ctx);
            }
        }
        /* Now there's an available spot in index i */
        vec->inds[i] = col;
        vec->nnz++;
        return gr_set(GR_ENTRY(vec->entries, i, sz), entry, ctx);
    }
    return GR_UNABLE; /* This shouldn't happen */
}
