
#include "acb_theta.h"

void
acb_theta_dupl_all_const(acb_ptr th2, acb_srcptr th, slong g, slong prec)
{
    acb_ptr v1, v2;
    acb_ptr res;
    slong n = 1<<g;
    ulong a, b;

    v1 = _acb_vec_init(n);
    v2 = _acb_vec_init(n);
    res = _acb_vec_init(n*n);

    for (a = 0; a < n; a++)
    {
        /* Set v2 to modified theta values */
        for (b = 0; b < n; b++)
        {
            acb_set(&v2[b], &th[b]);
            if (acb_theta_char_dot(a, b, g) == 1) acb_neg(&v2[b], &v2[b]);
        }
        acb_theta_agm_hadamard(v1, th, g, prec);
        acb_theta_agm_hadamard(v2, v2, g, prec);
        for (b = 0; b < n; b++) acb_mul(&v1[b], &v1[b], &v2[b], prec);
        acb_theta_agm_hadamard(v1, v1, g, prec);
        for (b = 0; b < n; b++) acb_mul_2exp_si(&res[n*a + b], &v1[b], -2*g);
    }    
    _acb_vec_set(th2, res, n*n);
    
    _acb_vec_clear(v1, n);
    _acb_vec_clear(v2, n);
    _acb_vec_clear(res, n*n);
}
