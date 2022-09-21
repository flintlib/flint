
#include "acb_theta.h"

void
acb_theta_agm_ctx_init_ext(acb_theta_agm_ctx_t ctx, acb_srcptr z,
        const acb_mat_t tau)
{
    slong g = acb_mat_nrows(tau);
    slong nb = 1<<g;
    slong dim = 2*(nb - 1);

    acb_theta_agm_ctx_is_ext(ctx) = 1;
    acb_theta_agm_ctx_dim(ctx) = dim;            
    acb_theta_agm_ctx_init_internal(ctx, nb, g);

    acb_mat_init(acb_theta_agm_ctx_tau(ctx), g, g);
    acb_theta_agm_ctx_z(ctx) = _acb_vec_init(2*g);
    acb_theta_agm_ctx_th(ctx) = _acb_vec_init(1<<(g+1));
    
    acb_mat_set(acb_theta_agm_ctx_tau(ctx), tau);
    _acb_vec_set(acb_theta_agm_ctx_z(ctx), z, g);
}
