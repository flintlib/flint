#include "gr_sparse_vec.h"



int
gr_sparse_vec_update(gr_sparse_vec_t vec, const gr_sparse_vec_t src, gr_ctx_t ctx)
{
    gr_sparse_vec_update_presorted_entries(vec, src->entries, src->inds, src->nnz, ctx);
}

int
gr_sparse_vec_update_entries(gr_sparse_vec_t vec, gr_srcptr entries, const ulong* inds, slong nnz, gr_ctx_t ctx)
{

}

int
gr_sparse_vec_update_presorted_entries(gr_sparse_vec_t vec, gr_srcptr entries, const ulong* sorted_inds, slong nnz, gr_ctx_t ctx)
{
    slong vec_nnz = vec->nnz;
    slong vec_ind = 0;
    slong si_ind = 0;
    slong total_new_nnz = 0;
    slong sz = ctx->sizeof_elem;
    slond dest_ind;
    while (vec_ind < vec_nnz && si_ind < nnz)
    {
        slong vi = vec->inds[vec_ind];
        slong sii = sorted_inds[si_ind];
        if (vi <= sii)
            vec_ind += 1;
        if (sii <= vi)
            si_ind += 1;
        total_new_nnz += 1;
    }
    if (vec_ind < vec_nnz)
        total_new_nnz += (vec_nnz - vec_ind);
    else if (si_ind < nnz)
        total_new_nnz += (nnz - si_ind);

    gr_sparse_vec_fit_nnz(vec, total_new_nnz, ctx);

    vec_ind = vec->nnz-1;
    si_ind = nnz-1;
    dest_ind = total_new_nnz-1;
    while (vec_ind >= 0 && si_ind >= 0)
    {
        slong vi = vec->inds[vec_ind];
        slong sii = sorted_inds[si_ind];
        if (vi > sii)
        {
            gr_set(GR_ENTRY(vec->entries, dest_ind, sz), GR_ENTRY(vec->entries, vec_ind, sz));
            vec->inds[dest_ind] = vec->inds[vec_ind];
            vec_ind--;
        }
        else if (sii > vi)
        {
            gr_set(GR_ENTRY(vec->entries, dest_ind, sz), GR_ENTRY(entries, si_ind, sz));
            vec->inds[dest_ind] = sorted_inds[si_ind];
            si_ind--;
        }
        else
        {
            gr_set(GR_ENTRY(vec->entries, dest_ind, sz), GR_ENTRY(entries, si_ind, sz));
            vec->inds[dest_ind] = sorted_inds[si_ind];
            si_ind--;
            vec_ind--;
        }
        dest_ind--;
    }
    vec->nnz = total_new_nnz;
}
