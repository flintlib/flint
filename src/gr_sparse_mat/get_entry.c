#include <stdlib.h>
#include "gr_sparse_mat.h"


static int ulong_cmp(const void* a, const void* b)
{
    ulong av = *((ulong*)(a));
    ulong bv = *((ulong*)(b));
    return (av < bv ? -1 : (av > bv ? 1 : 0));
}

gr_ptr
gr_csr_mat_find_entry(gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx)
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

int
gr_csr_mat_get_entry(gr_ptr res, gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    ulong* bs = NULL;

    if (row < 0 || row >= mat->r || col < 0 || col >= mat->c)
        return GR_DOMAIN;

    bs = bsearch(&col, mat->cols + mat->rows[row], mat->rows[row+1] - mat->rows[row], sizeof(ulong), ulong_cmp);

    if (bs == NULL)
        return gr_zero(res, ctx);
    return gr_set(res, GR_ENTRY(mat->nzs, bs - mat->cols, sz), ctx);
}

int
gr_lil_mat_get_entry(gr_ptr res, gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r)
        return GR_DOMAIN;

    return gr_sparse_vec_get_entry(res, &mat->rows[row], col, ctx);
}
