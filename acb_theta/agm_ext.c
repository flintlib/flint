
#include "acb_theta.h"

void acb_theta_agm_ext(acb_t r, acb_srcptr a, acb_srcptr all_r0, slong nb_bad,
		       slong nb_total, slong g, slong prec)
{  
  acb_ptr v;
  acb_t exp;  
  slong n = 1<<g;
  slong k;

  v = _acb_vec_init(2*n);
  acb_init(exp);
  
  _acb_vec_set(v, a, 2*n);
  
  for (k = 0; k < nb_bad; k++)
    {
      acb_theta_agm_ext_step_bad(v, v, all_r0 + k*2*n, g, prec);
    }
  for (k = nb_bad; k < nb_total; k++)
    {
      acb_theta_agm_ext_step_good(v, v, g, prec);
    }

  acb_div(r, &v[0], &v[n], prec);
  acb_one(exp);
  acb_mul_2exp_si(exp, exp, nb_total);
  acb_pow(r, r, exp, prec);

  _acb_vec_clear(v, n);
  acb_clear(exp);
}
