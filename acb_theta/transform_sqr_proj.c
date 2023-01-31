
#include "acb_theta.h"

void
acb_theta_transform_sqr_proj(acb_ptr res, acb_srcptr th2, const fmpz_mat_t mat,
                             slong prec)
{
    acb_ptr aux;
    slong g = fmpz_mat_nrows(mat) / 2;
    ulong n = 1 << g;
    ulong ab;
    ulong image_ab;
    fmpz_t eps;
    acb_t c;

    aux = _acb_vec_init(n);
    fmpz_init(eps);
    acb_init(c);

    for (ab = 0; ab < n; ab++)
    {
        image_ab = acb_theta_transform_image_char(eps, ab, mat);
        acb_unit_root(c, 4, prec);
        acb_pow_fmpz(c, c, eps, prec);
        acb_mul(c, c, &th2[image_ab], prec);
        acb_set(&aux[ab], c);
    }

    _acb_vec_set(res, aux, n);

    _acb_vec_clear(aux, n);
    fmpz_clear(eps);
    acb_clear(c);
}
