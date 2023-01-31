
#include "acb_theta.h"

void
fmpz_mat_get_b(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    slong j, k;

    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            fmpz_set(fmpz_mat_entry(res, j, k), fmpz_mat_entry(mat, j, k + g));
        }
    }
}
