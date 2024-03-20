#include <stdlib.h>
#include "gr_sparse_vec.h"

typedef struct
{
    slong i;
    slong ind;
}
sparse_vec_index_t;

static int sparse_vec_index_cmp(const void* a, const void* b)
{
    slong aind = ((sparse_vec_index_t*)(a))->ind;
    slong bind = ((sparse_vec_index_t*)(b))->ind;
    return (aind < bind ? -1 : (aind > bind ? 1 : 0));
}


static sparse_vec_index_t * _sort_inds(ulong * inds, slong num)
{
    slong i;
    sparse_vec_index_t * si;

    si = flint_malloc(num * sizeof(sparse_vec_index_t));
    for (i = 0; i < num; i++)
    {
        si[i].i = i;
        si[i].ind = inds[i];
    }

    qsort(si, num, sizeof(sparse_vec_index_t), sparse_vec_index_cmp);
    return si;
}

int
gr_sparse_vec_from_entries(gr_sparse_vec_t vec, ulong * inds, gr_srcptr entries, slong nnz, truth_t is_canonical, gr_ctx_t ctx)
{
    slong i;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    sparse_vec_index_t *si;
    gr_ptr vec_entry, entry;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);

    for (i = 0; i < nnz; ++i)
        if (inds[i] >= vec->length)
            return GR_DOMAIN;

    gr_sparse_vec_fit_nnz(vec, nnz, ctx);
    if (is_canonical == T_TRUE)
    {
        // Just copy data
        memcpy(vec->inds, inds, nnz * sizeof(ulong));
        status |= _gr_vec_set(vec->nzs, entries, nnz, ctx);
        vec->nnz = nnz;
    }
    else
    {
        si = _sort_inds(inds, nnz);
        vec->nnz = 0;
        vec_entry = NULL;
        for(i = 0; i < nnz; ++i)
        {
            entry = GR_ENTRY(entries, si[i].i, sz);
    
            // If index is repeated, accumulate into current vector entry
            if (i > 0 && si[i].ind == si[i-1].ind)
                status |= gr_add(vec_entry, vec_entry, entry, ctx);
            else
            {
                // If current entry is empty or nonzero, move to next one
                if (vec_entry == NULL || is_zero(vec_entry, ctx) != T_TRUE)
                {
                    vec_entry = GR_ENTRY(vec->nzs, vec->nnz, sz);
                    ++vec->nnz;
                }
                vec->inds[vec->nnz-1] = si[i].ind;
                status |= gr_set(vec_entry, entry, ctx);
            }
        }
        // Check if last entry accumulated to zero
        if (vec_entry != NULL && is_zero(vec_entry, ctx) == T_TRUE)
            vec->nnz--;

        flint_free(si);
    }
    return status;
}



