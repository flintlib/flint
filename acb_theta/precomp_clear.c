
#include "acb_theta.h"

void acb_theta_precomp_clear(acb_theta_precomp_t D)
{
  slong nb = D->nb;
  acb_mat_clear(acb_theta_precomp_exp_mat(D));
  flint_free(D->box);
  flint_free(D->indices);
  if (nb > 0) _acb_vec_clear(D->sqr_powers, nb);
}
