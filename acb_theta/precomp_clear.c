
#include "acb_theta.h"

void
acb_theta_precomp_clear(acb_theta_precomp_t D)
{
    slong g = acb_mat_nrows(acb_theta_precomp_exp_mat(D));
    slong nb_pow = D->indices[g];
    
    acb_mat_clear(acb_theta_precomp_exp_mat(D));
    flint_free(D->indices);
    if (nb_pow > 0) _acb_vec_clear(D->sqr_powers, nb_pow);
    _acb_vec_clear(D->exp_z, g * acb_theta_precomp_nb_z(D));
}
