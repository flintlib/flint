
#include "acb_theta.h"

void
acb_theta_dupl_all(acb_ptr th2, acb_srcptr th, slong g, slong prec)
{
    acb_ptr v1, v2, v3;
    acb_ptr res;
    slong n = 1 << g;
    ulong a, b;

    v1 = _acb_vec_init(n);
    v2 = _acb_vec_init(n);
    v3 = _acb_vec_init(n);
    res = _acb_vec_init(2 * n * n);

    for (a = 0; a < n; a++)
    {
        /* Set v3 to modified theta(0,.) values */
        for (b = 0; b < n; b++)
        {
            acb_set(&v3[b], &th[b + n]);
            if (acb_theta_char_dot(a, b, g) == 1)
                acb_neg(&v3[b], &v3[b]);
        }
        acb_theta_agm_hadamard(v1, th, g, prec);
        acb_theta_agm_hadamard(v2, th + n, g, prec);
        acb_theta_agm_hadamard(v3, v3, g, prec);
        for (b = 0; b < n; b++)
        {
            acb_mul(&v1[b], &v1[b], &v3[b], prec);
            acb_mul(&v2[b], &v2[b], &v3[b], prec);
        }
        acb_theta_agm_hadamard(v1, v1, g, prec);
        acb_theta_agm_hadamard(v2, v2, g, prec);
        for (b = 0; b < n; b++)
        {
            acb_mul_2exp_si(&res[n * a + b], &v1[b], -2 * g);
            acb_mul_2exp_si(&res[n * n + n * a + b], &v2[b], -2 * g);
        }
    }
    _acb_vec_set(th2, res, 2 * n * n);

    _acb_vec_clear(v1, n);
    _acb_vec_clear(v2, n);
    _acb_vec_clear(v3, n);
    _acb_vec_clear(res, 2 * n * n);
}
