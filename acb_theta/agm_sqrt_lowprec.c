
#include "acb_theta.h"

void
acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t x, const acb_t root, slong prec)
{
    acb_t res;
    acb_t dist;
  
    acb_init(res);
    acb_init(dist);

    /* Take any square root, avoiding potentially massive precision loss
       if x intersects the default branch cut */
  
    if (arb_contains_zero(acb_imagref(x))
	    && arb_is_negative(acb_realref(x)))
    {
	acb_neg(res, x);
	acb_sqrt(res, res, prec);
	acb_mul_onei(res, res);
    }
    else acb_sqrt(res, x, prec);

    /* Change sign if not contained in root */
    if (!acb_contains(root, res)) acb_neg(res, res);
    if (!acb_contains(root, res))
    {
	flint_printf("acb_theta_agm_sqrt_lowprec: Error (no suitable square root)\n");
	fflush(stdout);
	flint_abort();
    }
    
    acb_set(r, res);
    acb_clear(res);
    acb_clear(dist);
}
