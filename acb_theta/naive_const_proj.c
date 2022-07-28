
#include "acb_theta.h"

void acb_theta_naive_const_proj(acb_ptr th, const acb_mat_t tau, slong prec)
{
  slong g = acb_mat_nrows(tau);
  slong k;
  
  acb_theta_naive_const(th, tau, prec);
  for (k = 1; k < (1<<g); k++)
    {
      acb_div(&th[k], &th[k], &th[0], prec);
    }
  acb_one(&th[0]);
}
