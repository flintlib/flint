#include "gr_mat.h"

int
gr_mat_concat_vertical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;
    
    if (mat1->c != mat2->c || res->r != mat1->r + mat2->r)
        return GR_DOMAIN;

    if (c1 > 0)
    {
        for (i = 0; i < r1; i++)
            status |= _gr_vec_set(res->rows[i], mat1->rows[i], c1, ctx);
        for (i = 0; i < r2; i++)
            status |= _gr_vec_set(res->rows[i + r1], mat2->rows[i], c1, ctx);
    }

    return status;
}
