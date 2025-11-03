#include "fmpz.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_rand(fmpz_mod_mat_t mat, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    fmpz* e;
    slong i, j, r, c;
    r = fmpz_mod_mat_nrows(mat, ctx);
    c = fmpz_mod_mat_nrows(mat, ctx);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            e = fmpz_mod_mat_entry(mat, i, j);
            fmpz_randm(e, state, ctx->n);
        }
    }
}
