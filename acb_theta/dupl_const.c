
#include "acb_theta.h"

void
acb_theta_dupl_const(acb_ptr th2, acb_srcptr th, slong g, slong prec)
{
    acb_theta_agm_step_sqrt(th2, th, g, prec);
}
