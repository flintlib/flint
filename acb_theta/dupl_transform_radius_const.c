
#include "acb_theta.h"

void
acb_theta_dupl_transform_radius_const(arf_t rho, const arf_t r, acb_srcptr th,
                                      const fmpz_mat_t mat, slong prec)
{
    acb_ptr th_dupl;
    slong g = fmpz_mat_nrows(mat) / 2;
    slong n = 1 << g;

    th_dupl = _acb_vec_init(n * n);

    acb_theta_dupl_all_const(th_dupl, th, g, prec);
    acb_theta_transform_sqr_radius(rho, r, th_dupl, mat, prec);
    acb_theta_dupl_radius(rho, rho, th, n, prec);

    _acb_vec_clear(th_dupl, n * n);
}
