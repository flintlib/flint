
#include "acb_theta.h"

void acb_theta_naive_ind_const(acb_t th, ulong ab, const acb_mat_t tau, slong prec)
{
  slong g = acb_mat_nrows(tau);
  acb_ptr z;
  
  z = _acb_vec_init(g);

  _acb_vec_zero(z, g);
  acb_theta_naive_ind(th, ab, z, tau, prec);

  _acb_vec_clear(z, g);
}
