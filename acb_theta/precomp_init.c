
#include "acb_theta.h"

void
acb_theta_precomp_init(acb_theta_precomp_t D, slong g)
{
    acb_mat_init(acb_theta_precomp_exp_mat(D), g, g);
    D->indices = flint_malloc((g+1) * sizeof(slong));
    D->indices[g] = 0;
}
