
#include "acb_theta.h"

void
acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr z;
  
    z = _acb_vec_init(g);

    _acb_vec_zero(z, g);
    acb_theta_naive_all(th, z, 1, tau, prec);

    _acb_vec_clear(z, g);
}
