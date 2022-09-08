
#include "acb_theta.h"

void
acb_theta_precomp_clear(acb_theta_precomp_t D)
{
    slong g = acb_mat_nrows(acb_theta_precomp_exp_mat(D));
    slong nb = D->indices[g];
    
    acb_mat_clear(acb_theta_precomp_exp_mat(D));
    flint_free(D->indices);
    if (nb > 0) _acb_vec_clear(D->sqr_powers, nb);
}
