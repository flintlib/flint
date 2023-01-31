
#include "acb_theta.h"

void
acb_theta_agm_roots(acb_ptr roots, const acb_mat_t tau, slong nb_bad,
                    slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    slong k;

    acb_mat_init(w, g, g);
    acb_mat_set(w, tau);

    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_naive_const(roots + k * n, w, prec);
        acb_mat_scalar_mul_2exp_si(w, w, 1);
    }

    acb_mat_clear(w);
}
