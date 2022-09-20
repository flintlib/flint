
#include "acb_theta.h"

void
acb_theta_newton_const_half_proj(acb_ptr th, const acb_mat_t tau, slong prec)
{
    acb_theta_agm_ctx_t ctx;
    acb_mat_t half;
  
    slong g = acb_mat_nrows(tau);
    slong baseprec = ACB_THETA_AGM_BASEPREC;
    int stop = 0;
    int naive = 0;
  
    acb_theta_agm_ctx_init(ctx, tau);
    acb_mat_init(half, g, g);
    
    acb_mat_scalar_mul_2exp_si(half, tau, -1);

    /* Attempt to set up newton context */
    while (!stop)
    {
        acb_theta_agm_ctx_set(ctx, baseprec);
        if (!acb_theta_agm_ctx_is_valid(ctx))
	{
            baseprec *=2;
            if (baseprec > prec / ACB_THETA_AGM_BASEPREC_MAXQ)
	    {
                stop = 1;
                naive = 1;
	    }
	}
        else stop = 1;
    }

    if (naive) acb_theta_naive_const_proj(th, half, prec);
    else acb_theta_newton_run(th, ctx, prec);

    flint_printf("newton_const_half_proj:\n");
    for (stop = 0; stop < (1<<g); stop++)
    {
        acb_printd(&th[stop], 10); flint_printf("\n");
    }
    
    acb_mat_clear(half);
    acb_theta_agm_ctx_clear(ctx);
}
