#include <stdlib.h>
#include "gr_sparse_mat.h"

static int check_coords(ulong * rows, ulong * cols, slong num, slong r, slong c)
{
    slong i;

    for (i = 0; i < num; ++i)
        if (rows[i] >= r || cols[i] >= c)
            return 0;
    return 1;
}

int
gr_coo_mat_from_entries(gr_coo_mat_t mat, ulong *rows, ulong *cols, gr_srcptr entries, slong nnz, int is_canonical, gr_ctx_t ctx)
{
    int status;

    if (!check_coords(rows, cols, nnz, mat->r, mat->c))
        return GR_DOMAIN;

    mat->nnz = nnz;
    gr_coo_mat_fit_nnz(mat, nnz, ctx);
    memcpy(mat->rows, rows, nnz * sizeof(ulong));
    memcpy(mat->cols, cols, nnz * sizeof(ulong));
    status = _gr_vec_set(mat->nzs, entries, nnz, ctx);
    mat->is_canonical = is_canonical;
    return status;
}

