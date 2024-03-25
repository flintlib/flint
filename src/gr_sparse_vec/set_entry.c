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
        i = bs - vec->inds;
        if (gr_is_zero(entry, ctx) == T_TRUE)
        {
            // Shift everything above i down
            memmove(vec->inds + i, vec->inds + i + 1, (vec->nnz - i - 1)*sizeof(slong));
            for (j = i; j < vec->nnz; j++)
            {
                gr_swap(GR_ENTRY(vec->nzs, j, sz), GR_ENTRY(vec->nzs, j + 1, sz), ctx);
            }                        
            --vec->nnz;
            return GR_SUCCESS;
        }
    }
    else
    {
        if (gr_is_zero(entry, ctx) == T_TRUE)
        {
            // Already 0
            return GR_SUCCESS;
        }
        // Make room for new element
        gr_sparse_vec_fit_nnz(vec, vec->nnz+1, ctx);

        // Find location to put new element
        for (i = 0; i < nnz; i++)
            if (col < vec->inds[i])
                break;

        // Shift everything above i up
        memmove(vec->inds + i + 1, vec->inds + i, (vec->nnz - i)*sizeof(slong));
        for (j = vec->nnz; j > i; j--)
        {
            gr_swap(GR_ENTRY(vec->nzs, j-1, sz), GR_ENTRY(vec->nzs, j, sz), ctx);
        }
        vec->inds[i] = col;
        ++vec->nnz;
    }
    return gr_set(GR_ENTRY(vec->nzs, i, sz), entry, ctx);
}
