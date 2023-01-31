
#include "acb_theta.h"

void
acb_theta_agm_ctx_init_internal(acb_theta_agm_ctx_t ctx, slong nb, slong g)
{
    slong k;

    acb_mat_init(acb_theta_agm_ctx_tau(ctx), g, g);
    acb_theta_agm_ctx_z(ctx) = _acb_vec_init(2 * g);
    acb_theta_agm_ctx_th(ctx) = _acb_vec_init(1 << (g + 1));

    acb_theta_agm_ctx_nb(ctx) = nb;
    ctx->mat = flint_malloc(nb * sizeof(fmpz_mat_struct));
    ctx->k2 = flint_malloc(nb * sizeof(slong));
    ctx->ab = flint_malloc(nb * sizeof(slong));
    ctx->eps = flint_malloc(nb * sizeof(fmpz));
    ctx->nb_bad = flint_malloc(nb * sizeof(slong));
    for (k = 0; k < nb; k++)
        acb_theta_agm_ctx_nb_bad(ctx, k) = 0;
    ctx->roots = flint_malloc(nb * sizeof(acb_ptr));

    for (k = 0; k < nb; k++)
    {
        fmpz_mat_init(acb_theta_agm_ctx_mat(ctx, k), 2 * g, 2 * g);
        fmpz_init(acb_theta_agm_ctx_eps(ctx, k));
    }
}
