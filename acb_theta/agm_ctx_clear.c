
#include "acb_theta.h"

void acb_theta_agm_ctx_clear(acb_theta_agm_ctx_t ctx)
{
    slong k;
    slong n = acb_theta_agm_ctx_nb(ctx);

    acb_mat_clear(acb_theta_agm_ctx_tau(ctx));
    if (acb_theta_agm_ctx_is_ext(ctx))
    {
        _acb_vec_clear(acb_theta_agm_ctx_z(ctx), g);
        _acb_vec_clear(acb_theta_agm_ctx_th(ctx), 1<<(g+1));
    }
    else
    {
        _acb_vec_clear(acb_theta_agm_ctx_th(ctx), 1<<g);
    }        
  
    for (k = 0; k < n; k++) fmpz_mat_clear(acb_theta_agm_ctx_matrix(ctx, k));
    flint_free(ctx->matrices);
    flint_free(ctx->k2);
    flint_free(ctx->ab);
    for (k = 0; k < n; k++) fmpz_clear(acb_theta_agm_ctx_eps(ctx, k));
    flint_free(ctx->eps);
    for (k = 0; k < n; k++) acb_theta_agm_ctx_reset_steps(ctx, k, 0);
    flint_free(ctx->nb_bad_steps);
    flint_free(ctx->roots);
    flint_free(ctx->mi);  
    for (k = 0; k < n; k++) arf_clear(acb_theta_agm_ctx_M0(ctx, k));
    flint_free(ctx->M0);
    for (k = 0; k < n; k++) arf_clear(acb_theta_agm_ctx_minf(ctx, k));
    flint_free(ctx->minf);
    for (k = 0; k < n; k++) arf_clear(acb_theta_agm_ctx_rad(ctx, k));
    flint_free(ctx->rad);
    for (k = 0; k < n; k++) arf_clear(acb_theta_agm_ctx_max(ctx, k));
    flint_free(ctx->max);
    
    arf_clear(acb_theta_agm_ctx_rho(ctx));
    arf_clear(acb_theta_agm_ctx_M(ctx));
    arf_clear(acb_theta_agm_ctx_B3(ctx));
}
