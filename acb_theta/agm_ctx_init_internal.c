
#include "acb_theta.h"

void acb_theta_agm_ctx_init_internal(acb_theta_agm_ctx_t ctx, slong nb,
        slong g)
{
    slong k;
    
    acb_theta_agm_ctx_nb(ctx) = nb;    
    ctx->matrices = flint_malloc(nb * sizeof(fmpz_mat_struct));
    for (k = 0; k < nb; k++)
    {
        fmpz_mat_init(acb_theta_agm_ctx_matrix(ctx, k), 2*g, 2*g);
    }
    ctx->k2 = flint_malloc(nb * sizeof(slong));
    ctx->ab = flint_malloc(nb * sizeof(slong));
    ctx->eps = flint_malloc(nb * sizeof(fmpz));
    for (k = 0; k < nb; k++) fmpz_init(acb_theta_agm_ctx_eps(ctx, k));
    ctx->nb_bad_steps = flint_malloc(nb * sizeof(slong));
    for (k = 0; k < nb; k++) acb_theta_agm_ctx_nb_bad_steps(ctx, k) = 0;
    ctx->roots = flint_malloc(nb * sizeof(acb_ptr));
    ctx->mi = flint_malloc(nb * sizeof(arf_struct*));
    ctx->Mi = flint_malloc(nb * sizeof(arf_struct*));
    ctx->minf = flint_malloc(nb * sizeof(arf_struct));
    for (k = 0; k < nb; k++) arf_init(acb_theta_agm_ctx_minf(ctx, k));
    ctx->rad = flint_malloc(nb * sizeof(arf_struct));
    for (k = 0; k < nb; k++) arf_init(acb_theta_agm_ctx_rad(ctx, k));
    ctx->min = flint_malloc(nb * sizeof(arf_struct));
    for (k = 0; k < nb; k++) arf_init(acb_theta_agm_ctx_min(ctx, k));
    ctx->max = flint_malloc(nb * sizeof(arf_struct));
    for (k = 0; k < nb; k++) arf_init(acb_theta_agm_ctx_max(ctx, k));
    
    arf_init(acb_theta_agm_ctx_rho(ctx));
    arf_init(acb_theta_agm_ctx_M(ctx));
    arf_init(acb_theta_agm_ctx_B3(ctx));
}
