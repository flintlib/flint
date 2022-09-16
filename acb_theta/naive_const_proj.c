
#include "acb_theta.h"

void
acb_theta_naive_const_proj(acb_ptr th, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1<<g;
  
    acb_theta_naive_const(th, tau, prec);
    _acb_vec_scalar_div(&th[1], &th[1], n-1, &th[0], prec);
    acb_one(&th[0]);
}
