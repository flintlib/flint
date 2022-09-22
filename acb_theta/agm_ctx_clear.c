
#include "acb_theta.h"

void acb_theta_agm_ctx_clear(acb_theta_agm_ctx_t ctx)
{
    slong nb = acb_theta_agm_ctx_nb(ctx);
    slong g = acb_theta_agm_ctx_g(ctx);
    slong k;

    acb_mat_clear(acb_theta_agm_ctx_tau(ctx));
    _acb_vec_clear(acb_theta_agm_ctx_z(ctx), 2*g);
    _acb_vec_clear(acb_theta_agm_ctx_th(ctx), 1<<(g+1));

    for (k = 0; k < nb; k++)
    {
        fmpz_mat_clear(acb_theta_agm_ctx_matrix(ctx, k));
        fmpz_clear(acb_theta_agm_ctx_eps(ctx, k));
        acb_theta_agm_ctx_reset_steps(ctx, k, 0);
        arf_clear(acb_theta_agm_ctx_minf(ctx, k));
        arf_clear(acb_theta_agm_ctx_c(ctx, k));
        arf_clear(acb_theta_agm_ctx_c_ext(ctx, k));
        arf_clear(acb_theta_agm_ctx_e(ctx, k));
        arf_clear(acb_theta_agm_ctx_rad(ctx, k));
        arf_clear(acb_theta_agm_ctx_min(ctx, k));
        arf_clear(acb_theta_agm_ctx_max(ctx, k));
    }
    
    arf_clear(acb_theta_agm_ctx_rho(ctx));
    arf_clear(acb_theta_agm_ctx_M(ctx));
    arf_clear(acb_theta_agm_ctx_B3(ctx));
  
    flint_free(ctx->matrices);
    flint_free(ctx->k2);
    flint_free(ctx->ab);
    flint_free(ctx->eps);
    flint_free(ctx->nb_bad_steps);
    flint_free(ctx->roots);
    flint_free(ctx->mi);
    flint_free(ctx->Mi);
    flint_free(ctx->minf);
    flint_free(ctx->rad);
    flint_free(ctx->min);
    flint_free(ctx->max);
}
