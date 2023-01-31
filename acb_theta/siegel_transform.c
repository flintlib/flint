
#include "acb_theta.h"

void
acb_siegel_transform(acb_mat_t res, const fmpz_mat_t mat, const acb_mat_t tau,
                     slong prec)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    fmpz_mat_t a;
    acb_mat_t x, num, den;
    int r;

    fmpz_mat_init(a, g, g);
    acb_mat_init(x, g, g);
    acb_mat_init(num, g, g);
    acb_mat_init(den, g, g);

    fmpz_mat_get_a(a, mat);
    acb_mat_set_fmpz_mat(x, a);
    acb_mat_mul(num, x, tau, prec);
    fmpz_mat_get_b(a, mat);
    acb_mat_set_fmpz_mat(x, a);
    acb_mat_add(num, num, x, prec);

    acb_siegel_cocycle(den, mat, tau, prec);
    r = acb_mat_inv(den, den, prec);
    if (!r)
        acb_mat_indeterminate(den);

    acb_mat_mul(res, num, den, prec);

    fmpz_mat_clear(a);
    acb_mat_clear(x);
    acb_mat_clear(num);
    acb_mat_clear(den);
}
