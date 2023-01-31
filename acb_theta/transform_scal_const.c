
#include "acb_theta.h"

void
acb_theta_transform_scal_const(acb_t scal, const acb_mat_t tau,
                               const fmpz_mat_t mat, slong k2, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_t mu;
    acb_t det;
    acb_mat_t w;

    acb_init(mu);
    acb_init(det);
    acb_mat_init(w, g, g);

    acb_onei(mu);
    acb_pow_si(mu, mu, k2, prec);
    acb_siegel_cocycle(w, mat, tau, prec);
    acb_mat_det(det, w, prec);
    acb_mul(scal, det, mu, prec);

    acb_clear(mu);
    acb_clear(det);
    acb_mat_clear(w);
}
