
#include "acb_theta.h"

void
acb_theta_newton_eval(acb_ptr r, acb_srcptr th, const acb_theta_agm_ctx_t ctx,
        slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = acb_theta_agm_ctx_nb(ctx);
    acb_ptr dupl;
    acb_ptr transf;
    acb_ptr agm;
    arf_t err;
    slong nb_good;
    slong k, j;

    dupl = _acb_vec_init(1<<(2*g));
    transf = _acb_vec_init(1<<g);
    agm = _acb_vec_init(n);
    arf_init(err);
  
    acb_theta_dupl_all_const(dupl, th, g, prec);
    for (k = 0; k < n; k++)
    {      
        acb_theta_transform_sqr_proj(transf, dupl,
                acb_theta_agm_ctx_matrix(ctx, k), prec);
        for (j = 1; j < n; j++)
        {
            acb_div(&transf[j], &transf[j], &transf[0], prec);
        }
        acb_one(&transf[0]);
        nb_good = acb_theta_agm_nb_good_steps(err, g, prec);
        acb_theta_agm(&agm[k], transf, acb_theta_agm_ctx_roots(ctx, k), err,
                acb_theta_agm_ctx_nb_bad_steps(ctx, k), nb_good, g, prec);
    }
    for (k = 0; k < n-1; k++)
    {
        acb_div(&r[k], &agm[0], &agm[k+1], prec);
    }

    _acb_vec_clear(dupl, 1<<(2*g));
    _acb_vec_clear(transf, 1<<g);
    _acb_vec_clear(agm, n);
    arf_clear(err);
}
