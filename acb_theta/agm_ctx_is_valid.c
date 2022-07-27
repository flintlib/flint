
#include "acb_theta.h"

int acb_theta_agm_ctx_is_valid(const acb_theta_agm_ctx_t ctx)
{
  return arf_is_positive(acb_theta_newton_rho(ctx))
    && arf_is_finite(acb_theta_newton_max(ctx))
    && arf_is_finite(acb_theta_newton_inv_der(ctx));
}
