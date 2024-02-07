#include "gr_sparse_vec.h"

// int
// gr_vec_set(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx)
// {
//     int status = GR_SUCCESS;
//     if (res != src)
//     {
//         slong nnz = src->nnz;
//         gr_sparse_vec_fit_nnz(res, nnz, ctz);
//         memcpy(res->inds, src->inds, nnz*sizeof(slong));
//         res->nnz = nnz;
//         status = _gr_vec_set(res->entries, src->entries, nnz, ctx);
//     }
//     return status;
// }

int
gr_sparse_vec_update(gr_sparse_vec_t vec, const gr_sparse_vec_t src, gr_ctx_t ctx)
{
    gr_sparse_vec_update_entries(vec, src->entries, src->inds, src->nnz, ctx);
}

int
gr_sparse_vec_update_entries(gr_sparse_vec_t vec, gr_srcptr entries, const ulong* inds, slong nnz, gr_ctx_t ctx)
{

}