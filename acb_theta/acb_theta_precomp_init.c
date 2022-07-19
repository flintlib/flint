
#include "acb_theta.h"

void acb_theta_precomp_init(acb_theta_precomp_t D, slong g)
{
  acb_theta_precomp_g(D) = g;
  acb_mat_init(acb_theta_precomp_exp_mat(D), g, g);
  D->box = flint_malloc(g * sizeof(slong));
  D->indices = flint_malloc(g * sizeof(slong));
  D->nb = 0;
}
