#include "fmpz.h"
#include "fmpz_mod_mat.h"

/* entries uniformly random in [1, n), or 0 if n = 1 */
void fmpz_mod_mat_randfull(fmpz_mod_mat_t mat, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    fmpz* e;
    fmpz_t m;
    slong i, j, r, c;
    r = fmpz_mod_mat_nrows(mat, ctx);
    c = fmpz_mod_mat_ncols(mat, ctx);

    fmpz_init(m);
    fmpz_sub_ui(m, ctx->n, 1);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            e = fmpz_mod_mat_entry(mat, i, j);
            if (fmpz_is_zero(m))
            {
                fmpz_zero(e);
            }
            else
            {
                fmpz_randm(e, state, m);
                fmpz_add_ui(e, e, 1);
            }
        }
    }

    fmpz_clear(m);
}
