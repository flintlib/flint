#include "gr_sparse_mat.h"

void
gr_csr_mat_set_cols(gr_csr_mat_t mat, slong c, gr_ctx_t ctx)
{
    slong idx, new_idx, row;
    slong sz = ctx->sizeof_elem;

    if (c < mat->c)
    {
        // Keep track of row to update offsets
        row = 1;

        // May have some columns to discard
        for (idx = new_idx = 0; idx < mat->nnz; ++idx)
        {
            if (mat->rows[row] == idx)
            {
                mat->rows[row++] = new_idx;
            }
            if (mat->cols[idx] < c && idx != new_idx)
            {
                // Move this entry down in array
                mat->cols[new_idx] = mat->cols[idx];
                gr_swap(GR_ENTRY(mat->nzs, new_idx, sz), GR_ENTRY(mat->nzs, idx, sz), ctx);
                ++new_idx;
            }
        }
        mat->nnz = new_idx;
    }
   mat->c = c;
 }

void
gr_lil_mat_set_cols(gr_lil_mat_t mat, slong c, gr_ctx_t ctx)
{
    slong row;

    mat->c = c;
    for (row = 0; row < mat->r; ++row)
    {
        gr_sparse_vec_set_length(&mat->rows[row], c, ctx);
    }
}

void
gr_coo_mat_set_cols(gr_coo_mat_t mat, slong c, gr_ctx_t ctx)
{
    slong idx, new_idx;
    slong sz = ctx->sizeof_elem;

    if (c < mat->c)
    {
        // May have some columns to discard
        for (idx = new_idx = 0; idx < mat->nnz; ++idx)
        {
            if (mat->cols[idx] < c && idx != new_idx)
            {
                // Move this entry down in array
                mat->rows[new_idx] = mat->rows[idx];
                mat->cols[new_idx] = mat->cols[idx];
                gr_swap(GR_ENTRY(mat->nzs, new_idx, sz), GR_ENTRY(mat->nzs, idx, sz), ctx);
                ++new_idx;
            }
        }
        mat->nnz = new_idx;
    }
    mat->c = c;
 }
