#include "gr_sparse_vec.h"

int
gr_vec_set_sparse_vec(gr_ptr vec, gr_sparse_vec_t src, gr_ctx_t ctx)
{
    slong i,j,sz,nnz,status,len;
    status = GR_SUCCESS;
    sz = ctx->sizeof_elem;
    len = src->length;
    nnz = src->nnz;
    j = 0;
    for (i = 0; i < len; i++)
    {
        /* if i is the column of the next nonzero in src, copy it in */
        if (j < nnz && i == src->inds[j])
        {
            status |= gr_set(GR_ENTRY(vec, i, sz), GR_ENTRY(src->nzs, j, sz), ctx);
            j++;
        }
        else
        {
            status |= gr_zero(GR_ENTRY(vec, i, sz), ctx);
        }
    }
    return status;
}
