#include "gr_sparse_vec.h"

void
gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx)
{
    vec->length = len;
    /* Scan backward through the nonzeros to discard the ones off the end */
    /* Note that we don't actually free anything; we just mark the nonzeros as unused */
    slong i=vec->nnz-1;
    while (i>=0 && vec->cols[i] >= len) 
    {
        i--;
    }
    vec->nnz = i+1;
}
