
#include "acb_theta.h"

void acb_theta_vecsqr(acb_ptr th2, acb_srcptr th, slong n, slong prec)
{
    slong k;
    for (k = 0; k < n; k++)
    {
	acb_sqr(&th2[k], &th[k], prec);
    }
}
