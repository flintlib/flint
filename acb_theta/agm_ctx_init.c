
#include "acb_theta.h"

void acb_theta_agm_ctx_init(acb_theta_agm_ctx_t ctx, slong g, slong n)
{
  slong k;

  acb_theta_agm_ctx_g(ctx) = g;
  acb_theta_agm_ctx_nb(ctx) = n;
  ctx->matrices = flint_malloc(n * sizeof(fmpz_mat_struct));
  for (k = 0; k < n; k++) fmpz_mat_init(acb_theta_agm_ctx_matrix(ctx, k), 2*g, 2*g);
  ctx->nb_bad_steps = flint_malloc(n * sizeof(slong));
  for (k = 0; k < n; k++) acb_theta_agm_ctx_nb_bad_steps(ctx, k) = 0;
  ctx->roots = flint_malloc(n * sizeof(acb_ptr));
  ctx->mi = flint_malloc(n * sizeof(arf_struct*));
  ctx->M0 = flint_malloc(n * sizeof(arf_struct));
  for (k = 0; k < n; k++) arf_init(acb_theta_agm_ctx_M0(ctx, k));
  ctx->minf = flint_malloc(n * sizeof(arf_struct));
  for (k = 0; k < n; k++) arf_init(acb_theta_agm_ctx_minf(ctx, k));
  arf_init(acb_theta_agm_ctx_rho(ctx));
  arf_init(acb_theta_agm_ctx_max(ctx));
  arf_init(acb_theta_agm_ctx_inv_der(ctx));
}
