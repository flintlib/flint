
#include "acb_theta.h"

int
arb_mat_is_nonsymmetric(const arb_mat_t mat)
{
    arb_mat_t tp;
    slong nrows = arb_mat_nrows(mat);
    int res;

    if (nrows != arb_mat_ncols(mat))
        return 1;

    arb_mat_init(tp, nrows, nrows);
    arb_mat_transpose(tp, mat);
    res = !arb_mat_overlaps(tp, mat);
    arb_mat_clear(tp);

    return res;
}
