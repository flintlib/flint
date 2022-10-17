
#include "acb_theta.h"

int
acb_theta_is_balanced(slong* j0, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_t test;
    slong j;
    int r = 1;

    arb_init(test);
    
    for (j = 0; j < g-1; j++)
    {
	arb_mul_si(test, acb_imagref(arb_mat_entry(tau, j, j)),
		ACB_THETA_BALANCE_THRESHOLD, prec);
	if (arb_lt(test, acb_imagref(arb_mat_entry(tau, j+1, j+1))))
	{	    
	    r = 0;
	    *j0 = j;
	    break;
	}
    }

    arb_clear(test);
    return r;
}
