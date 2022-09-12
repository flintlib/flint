
#include "acb_theta.h"

int
acb_theta_agm_ctx_is_valid(const acb_theta_agm_ctx_t ctx)
{
    return arf_cmp_si(acb_theta_agm_ctx_rho(ctx), 0) > 0
        && arf_is_finite(acb_theta_agm_ctx_max(ctx))
        && arf_is_finite(acb_theta_agm_ctx_inv_der(ctx));
}
