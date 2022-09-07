
#include "acb_theta.h"

void
acb_theta_agm_step_sqrt(acb_ptr r, acb_srcptr a, slong g, slong prec)
{
    acb_ptr v;
    slong n = 1<<g;
    slong k;

    v = _acb_vec_init(n);
  
    acb_theta_agm_hadamard(v, a, g, prec);  
    for (k = 0; k < n; k++) acb_sqr(&v[k], &v[k], prec);
    acb_theta_agm_hadamard(r, v, g, prec);
    for (k = 0; k < n; k++) acb_mul_2exp_si(&r[k], &r[k], -2*g);

    _acb_vec_clear(v, n);
}
