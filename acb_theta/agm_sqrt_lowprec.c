
#include "acb_theta.h"

void
acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t a, const acb_t root,
        slong prec)
{
    acb_t res;
    acb_t dist;
  
    acb_init(res);
    acb_init(dist);

    /* Take any square root, avoiding potentially massive precision loss
       if a intersects the default branch cut */
  
    if (arb_contains_zero(acb_imagref(a))
	    && arb_is_negative(acb_realref(a)))
    {
	acb_neg(res, a);
	acb_sqrt(res, res, prec);
	acb_mul_onei(res, res);
    }
    else acb_sqrt(res, a, prec);

    /* Change sign if not contained in root */
    if (!acb_overlaps(root, res)) acb_neg(res, res);
    if (!acb_overlaps(root, res))
    {
	flint_printf("acb_theta_agm_sqrt_lowprec: Error (no suitable square root)\n");
        acb_printd(root, 10); flint_printf("\n");
        acb_printd(res, 10); flint_printf("\n");
	fflush(stdout);
	flint_abort();
    }
    
    acb_set(r, res);
    acb_clear(res);
    acb_clear(dist);
}
