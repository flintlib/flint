
#include "acb_theta.h"

void
acb_siegel_transform_ext(acb_ptr r, acb_mat_t w, const fmpz_mat_t mat,
                         acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    fmpz_mat_t a;
    acb_mat_t x, num, den;
    acb_mat_t vec;
    int res;
    slong k;

    fmpz_mat_init(a, g, g);
    acb_mat_init(x, g, g);
    acb_mat_init(num, g, g);
    acb_mat_init(den, g, g);
    acb_mat_init(vec, g, 1);

    fmpz_mat_get_a(a, mat);
    acb_mat_set_fmpz_mat(x, a);
    acb_mat_mul(num, x, tau, prec);
    fmpz_mat_get_b(a, mat);
    acb_mat_set_fmpz_mat(x, a);
    acb_mat_add(num, num, x, prec);

    acb_siegel_cocycle(den, mat, tau, prec);
    res = acb_mat_inv(den, den, prec);
    if (!res)
        acb_mat_indeterminate(den);

    acb_mat_mul(w, num, den, prec);

    acb_mat_transpose(den, den);
    for (k = 0; k < g; k++)
        acb_set(acb_mat_entry(vec, k, 0), &z[k]);
    acb_mat_mul(vec, den, vec, prec);
    for (k = 0; k < g; k++)
        acb_set(&r[k], acb_mat_entry(vec, k, 0));

    fmpz_mat_clear(a);
    acb_mat_clear(x);
    acb_mat_clear(num);
    acb_mat_clear(den);
    acb_mat_clear(vec);
}
