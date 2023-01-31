
#include "acb_theta.h"

void
acb_theta_newton_half_proj(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
                           slong prec)
{
    acb_theta_agm_ctx_t ctx;
    acb_mat_t half;
    acb_ptr z_0;

    slong g = acb_mat_nrows(tau);
    slong baseprec = ACB_THETA_AGM_BASEPREC;
    int stop = 0;
    int naive = 0;

    acb_theta_agm_ctx_init_ext(ctx, z, tau);
    acb_mat_init(half, g, g);
    z_0 = _acb_vec_init(2 * g);

    acb_mat_scalar_mul_2exp_si(half, tau, -1);

    /* Attempt to set up newton context */
    while (!stop)
    {
        stop = acb_theta_agm_ctx_set(ctx, baseprec);
        if (!stop)
        {
            baseprec *= 2;
            if (baseprec > prec / ACB_THETA_AGM_BASEPREC_MAXQ)
            {
                stop = 1;
                naive = 1;
            }
        }
    }

    if (naive)
    {
        _acb_vec_set(z_0, z, g);
        acb_theta_naive_proj(th, z_0, 2, half, prec);
    }
    else
    {
        acb_theta_newton_run(th, ctx, prec);
    }

    acb_theta_agm_ctx_clear(ctx);
    acb_mat_clear(half);
    _acb_vec_clear(z_0, 2 * g);
}
