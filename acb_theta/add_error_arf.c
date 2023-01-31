
#include "acb_theta.h"

void
arb_mat_add_error_arf(arb_mat_t mat, const arf_t err)
{
    slong k = acb_mat_nrows(mat);
    slong n = acb_mat_ncols(mat);
    slong i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < n; j++)
        {
            arb_add_error_arf(arb_mat_entry(mat, i, j), err);
        }
    }
}
