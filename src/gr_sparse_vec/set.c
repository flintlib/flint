#include "gr_sparse_vec.h"

int
gr_sparse_vec_set(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    if (dst != src)
    {
        slong nnz = src->nnz;
        gr_sparse_vec_fit_nnz(dst, nnz, ctx);
        memcpy(dst->inds, src->inds, nnz*sizeof(slong));
        status = _gr_vec_set(dst->nzs, src->nzs, nnz, ctx);
        dst->nnz = nnz;
        dst->length = src->length;
    }
    return status;
}
