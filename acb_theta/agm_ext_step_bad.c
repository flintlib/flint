
#include "acb_theta.h"

void acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr r0, slong g, slong prec)
{
  slong k;  
  for (k = 0; k < (1<<(g+1)); k++)
    {
      acb_theta_agm_sqrt_lowprec(&r[k], &a[k], &r0[k], prec);
    }
  acb_theta_agm_ext_step_sqrt(r, r, g, prec);
}
