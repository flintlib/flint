
#include "acb_theta.h"

void acb_theta_agm(acb_t r, acb_srcptr a, acb_srcptr all_r0, slong nb_bad,
		   slong nb_total, slong g, slong prec)
{
  acb_ptr v;
  slong n = 1<<g;
  slong k;

  v = _acb_vec_init(n);
  _acb_vec_set(v, a, n);
  
  for (k = 0; k < nb_bad; k++)
    {
      acb_theta_agm_step_bad(v, v, all_r0 + k*n, g, prec);
    }
  for (k = nb_bad; k < nb_total; k++)
    {
      if (k % ACB_THETA_CHECK_CV == 0)
	{
	  if (acb_theta_agm_is_reached(v, prec)) break;
	}
      acb_theta_agm_step_good(v, v, g, prec);
    }

  acb_set(r, &v[0]);
  _acb_vec_clear(v, n);
}
