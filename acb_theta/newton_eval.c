
#include "acb_theta.h"

void
acb_theta_newton_eval(acb_ptr r, acb_srcptr th,
        const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = acb_theta_agm_ctx_nb(ctx);
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);
    acb_ptr dupl;
    acb_ptr transf;
    acb_ptr agm;
    acb_t scal, scal2;
    arf_t err;
    slong k;

    if (is_ext)
    {        
        dupl = _acb_vec_init(1<<(2*g+1));
        transf = _acb_vec_init(1<<(g+1));
        agm = _acb_vec_init(2*n);
    }
    else
    {
        dupl = _acb_vec_init(1<<(2*g));
        transf = _acb_vec_init(1<<g);
        agm = _acb_vec_init(n);
    }
    acb_init(scal);
    acb_init(scal2);
    arf_init(err);

    /* Duplicate */
    if (is_ext)
    {
        acb_theta_dupl_all(dupl, th, g, prec);
    }
    else
    {
        acb_theta_dupl_all_const(dupl, th, g, prec);
    }
        
    /* Compute agms for each matrix */
    for (k = 0; k < n; k++)
    {
        /* Transform theta values */
        acb_theta_transform_sqr_proj(transf, dupl,
            acb_theta_agm_ctx_mat(ctx, k), prec);
        if (is_ext)
        {
            acb_theta_transform_sqr_proj(&transf[1<<g], &dupl[1<<(2*g)],
                    acb_theta_agm_ctx_mat(ctx, k), prec);
        }

        /* Projectivize */
        acb_set(scal, &transf[0]);
        _acb_vec_scalar_div(&transf[1], &transf[1], (1<<g)-1, scal, prec);
        acb_one(&transf[0]);
        if (is_ext)
        {   
            acb_set(scal, &transf[1<<g]);
            _acb_vec_scalar_div(&transf[(1<<g)+1], &transf[(1<<g)+1], (1<<g)-1,
                    scal, prec);
            acb_one(&transf[1<<g]);
        }

        /* Get agm */
        if (is_ext)
        {            
            acb_theta_agm_ext(&agm[k], &agm[n+k], transf,
                    acb_theta_agm_ctx_roots(ctx, k),
                    acb_theta_agm_ctx_nb_bad(ctx, k), g, prec);
        }
        else
        {
            acb_theta_agm(&agm[k], transf, acb_theta_agm_ctx_roots(ctx, k),
                    acb_theta_agm_ctx_nb_bad(ctx, k), g, prec);
        }
    }
    
    /* Renormalize dupl, as first matrix is I */
    acb_mul(scal, &agm[0], &dupl[0], prec);
    acb_inv(scal, scal, prec);
    if (is_ext)
    {        
        acb_mul(scal2, &agm[n], &dupl[1<<(2*g)], prec);
        acb_inv(scal2, scal2, prec);
    }
    
    for (k = 0; k < n-1; k++)
    {
        acb_mul(&r[k], &agm[k+1], &dupl[acb_theta_agm_ctx_ab(ctx, k+1)], prec);
        acb_mul(&r[k], scal, &r[k], prec);
        if (is_ext)
        {            
            acb_mul(&r[k+n-1], &agm[k+n+1],
                    &dupl[acb_theta_agm_ctx_ab(ctx, k+1) + (1<<(2*g))], prec);
            acb_mul(&r[k+n-1], scal2, &r[k+n-1], prec);
        }
    }

    flint_printf("newton_eval result:\n");    
    for (k = 0; k < acb_theta_agm_ctx_dim(ctx); k++)
    {
        acb_printd(&r[k], 10); flint_printf("\n");
    }

    /* Clear */
    if (is_ext)
    {        
        _acb_vec_clear(dupl, 1<<(2*g+1));
        _acb_vec_clear(transf, 1<<(g+1));
        _acb_vec_clear(agm, 2*n);
    }
    else
    {
        _acb_vec_clear(dupl, 1<<(2*g));
        _acb_vec_clear(transf, 1<<g);
        _acb_vec_clear(agm, n);
    }
    acb_clear(scal);
    acb_clear(scal2);
    arf_clear(err);
}
