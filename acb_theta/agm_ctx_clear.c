
#include "acb_theta.h"

void acb_theta_agm_ctx_clear(acb_theta_agm_ctx_t ctx)
{
  slong k, j;
  slong n = acb_theta_agm_ctx_nb(ctx);
  slong m;
  
  for (k = 0; k < n; k++) fmpz_mat_clear(acb_theta_agm_ctx_matrix(ctx, k));
  flint_free(ctx->matrices);
  for (k = 0; k < n; k++) acb_theta_agm_ctx_reset_steps(ctx, k, 0);
  flint_free(ctx->nb_bad_steps);
  flint_free(ctx->roots);
  flint_free(ctx->mi);  
  for (k = 0; k < n; k++) arf_clear(acb_theta_agm_ctx_M0(ctx, k));
  flint_free(ctx->M0);
  for (k = 0; k < n; k++) arf_clear(acb_theta_agm_ctx_minf(ctx, k));
  flint_free(ctx->minf);
  arf_clear(acb_theta_agm_ctx_rho(ctx));
  arf_clear(acb_theta_agm_ctx_max(ctx));
  arf_clear(acb_theta_agm_ctx_inv_der(ctx));
}
