#include <stdlib.h>
#include "gr_sparse_mat.h"


static int ulong_cmp(const void* a, const void* b)
{
    ulong av = *((ulong*)(a));
    ulong bv = *((ulong*)(b));
    return (av < bv ? -1 : (av > bv ? 1 : 0));
}

gr_ptr gr_csr_mat_find_entry(gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    ulong* bs = NULL;

    if (row < 0 || row >= mat->r || col < 0 || col >= mat->c)
        return NULL;

    bs = bsearch(&col, mat->cols + mat->rows[row], mat->rows[row+1] - mat->rows[row], sizeof(ulong), ulong_cmp);

    if (bs == NULL)
        return NULL;
    return GR_ENTRY(mat->nzs, bs - mat->cols, sz);
}

gr_ptr
gr_lil_mat_find_entry(gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r)
        return NULL;

    return gr_sparse_vec_find_entry(&mat->rows[row], col, ctx);
}

gr_ptr gr_coo_mat_find_entry(gr_coo_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    slong i, lower, upper;
    slong sz = ctx->sizeof_elem;

    if (row < 0 || row >= mat->r || col < 0 || col >= mat->c)
        return NULL;

    if (mat->is_canonical)
    {
        // Find using (manual) binary search
        lower = 0;
        upper = mat->nnz;
        while (lower <= upper)
        {
            i = (lower + upper)/2;
            if (mat->rows[i] < row || (mat->rows[i] == row && mat->cols[i] < col))
                lower = i + 1;
            else if (mat->rows[i] > row || (mat->rows[i] == row && mat->cols[i] > col))
                upper = i - 1;
            else
                return GR_ENTRY(mat->nzs, i, sz);
        }
    }
    else
    {
        // Exhaust to find first instance
        for (i = 0; i < mat->nnz; ++i)
            if (mat->rows[i] == row && mat->cols[i] == col)
                return GR_ENTRY(mat->nzs, i, sz);
    }
    return NULL;
}

int gr_csr_mat_get_entry(gr_ptr res, gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    gr_ptr ptr = gr_csr_mat_find_entry(mat, row, col, ctx);

    if (ptr == NULL)
        return gr_zero(res, ctx);
    return gr_set(res, ptr, ctx);
}

int gr_lil_mat_get_entry(gr_ptr res, gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    gr_ptr ptr = gr_lil_mat_find_entry(mat, row, col, ctx);

    if (ptr == NULL)
        return gr_zero(res, ctx);
    return gr_set(res, ptr, ctx);
}

int gr_coo_mat_get_entry(gr_ptr res, gr_coo_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    gr_ptr ptr = gr_coo_mat_find_entry(mat, row, col, ctx);

    if (ptr == NULL)
        return gr_zero(res, ctx);
    return gr_set(res, ptr, ctx);
}
