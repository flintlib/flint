
#include "acb_theta.h"

void
acb_theta_dupl_transform_radius(arf_t rho, const arf_t r, acb_srcptr th,
        const fmpz_mat_t mat, slong prec)
{
    slong g = fmpz_mat_nrows(mat)/2;
    acb_ptr th_dupl;
    slong n = 1<<g;
    arf_t aux;

    th_dupl = _acb_vec_init(2*n*n);
    arf_init(aux);
    
    acb_theta_dupl_all(th_dupl, th, g, prec);
    acb_theta_transform_sqr_radius(aux, r, th_dupl, mat, prec);
    acb_theta_transform_sqr_radius(rho, r, th_dupl + n*n, mat, prec);
    arf_min(rho, rho, aux);
    acb_theta_dupl_radius(rho, rho, th, 2*n, prec);
    
    _acb_vec_clear(th_dupl, 2*n*n);
    arf_clear(aux);
}
