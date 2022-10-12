
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
        fmpz_mat_clear(acb_theta_agm_ctx_mat(ctx, k));
        fmpz_clear(acb_theta_agm_ctx_eps(ctx, k));
        acb_theta_agm_ctx_reset_steps(ctx, k, 0);
    }
      
    flint_free(ctx->mat);
    flint_free(ctx->k2);
    flint_free(ctx->ab);
    flint_free(ctx->eps);
    flint_free(ctx->nb_bad);
    flint_free(ctx->roots);
}
