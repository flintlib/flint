
#include "acb_theta.h"

void acb_theta_agm_ext_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec)
{
  acb_ptr v;
  slong n = 1<<g;
  slong k;

  v = _acb_vec_init(2*n);
  
  acb_theta_agm_hadamard(v+n, a+n, g, prec);
  acb_theta_agm_hadamard(v, a, g, prec);
  for (k = 0; k < n; k++)
    {
      acb_mul(&v[k], &v[k], &v[k+n], prec);
      acb_sqr(&v[k+n], &v[k+n], prec);
    }
  
  acb_theta_agm_hadamard(r, v, g, prec);
  acb_theta_agm_hadamard(r, v+n, g, prec);
  for (k = 0; k < 2*n; k++)
    {
      acb_mul_2exp_si(&r[k], &r[k], -g);
    }

  _acb_vec_clear(v, 2*n);
}
