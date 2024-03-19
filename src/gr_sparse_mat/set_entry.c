#include <stdlib.h>
#include "gr_sparse_mat.h"

static int ulong_cmp(const void* a, const void* b)
{
    ulong av = *((ulong*)(a));
    ulong bv = *((ulong*)(b));
    return (av < bv ? -1 : (av > bv ? 1 : 0));
}

int
gr_csr_mat_set_entry(gr_csr_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)
{
    slong i, j;
    slong sz = ctx->sizeof_elem;
    ulong* bs = NULL;

    if (row < 0 || row >= mat->r || col < 0 || col >= mat->c)
        return GR_DOMAIN;

    bs = bsearch(&col, mat->cols + mat->rows[row], mat->rows[row+1] - mat->rows[row], sizeof(ulong), ulong_cmp);

    if (bs != NULL)
    {
        i = bs - mat->cols;
        if (gr_is_zero(entry, ctx) == T_TRUE)
        {
            // Shift everything above i down
            for (j = row; j <= mat->r; ++j)
                --mat->rows[j];
            memmove(mat->cols + i, mat->cols + i + 1, (mat->nnz - i - 1)*sizeof(slong));
            for (j = i; j < mat->nnz; j++)
            {
                gr_swap(GR_ENTRY(mat->nzs, j, sz), GR_ENTRY(mat->nzs, j + 1, sz), ctx);
            }                        
            --mat->nnz;

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
        // Make space for one more nonzero
        gr_csr_mat_fit_nnz(mat, mat->nnz+1, ctx);

        // Find location
        for (i = mat->rows[row]; i < mat->rows[row + 1]; ++i)
            if (col < mat->cols[i])
                break;

        // Shift everything above i up
        for (j = row; j <= mat->r; ++j)
            ++mat->rows[j];
        memmove(mat->cols + i + 1, mat->cols + i, (mat->nnz - i)*sizeof(slong));
        for (j = mat->nnz; j > i; j--)
        {
            gr_swap(GR_ENTRY(mat->nzs, j-1, sz), GR_ENTRY(mat->nzs, j, sz), ctx);
        }
        
        mat->cols[i] = col;
        ++mat->nnz;
    }
    return gr_set(GR_ENTRY(mat->nzs, i, sz), entry, ctx);
}

int
gr_lil_mat_set_entry(gr_lil_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r)
        return GR_DOMAIN;

    return gr_sparse_vec_set_entry(&mat->rows[row], col, entry, ctx);
}
