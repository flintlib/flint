
#include "acb_theta.h"

void
acb_theta_agm_ctx_update_bounds(acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong nb_th;
    slong nb_bad;
    arb_t abs;
    arb_t m;
    slong i, j;

    arb_init(abs);
    arb_init(m);

    nb_th = n;
    if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;

    nb_bad = acb_theta_agm_ctx_nb_bad_steps(ctx, k);
    
    /* Compute M0 */
    arb_zero(m);
    for (i = 0; i < nb_th * nb_bad; i++)
    {
        acb_abs(abs, &acb_theta_agm_ctx_roots(ctx, k)[i], prec);
        arb_max(m, m, abs, prec);
    }
    arb_sqr(m, m, prec);
    arb_get_ubound_arf(acb_theta_agm_ctx_M0(ctx, k), m, prec);
    
    /* Compute minf */
    acb_abs(m, &acb_theta_agm_ctx_roots(ctx, k)[nb_th * (nb_bad-1)], prec);
    arb_div_si(m, m, 20, prec); /* Cf. agm_nb_bad_steps */
    arb_get_lbound_arf(acb_theta_agm_ctx_minf(ctx, k), m, prec);
    
    /* Compute mi */
    for (i = 0; i < nb_bad; i++)
    {
        arb_pos_inf(m);
        for (j = 0; j < nb_th; j++)
        {
            acb_abs(abs, &acb_theta_agm_ctx_roots(ctx,k)[i*nb_th + j], prec);
            arb_min(m, m, abs, prec);
        }
        arb_sqr(m, m, prec);
        arb_get_lbound_arf(&acb_theta_agm_ctx_mi(ctx,k)[i], m, prec);
    }

    arb_clear(abs);
    arb_clear(m);
}
