
#include "acb_theta.h"

/* Collect data for a given symplectic matrix */

static void
set_matrix(acb_theta_agm_ctx_t ctx, slong k, const acb_mat_t tau,
        const fmpz_mat_t N, slong prec)
{
    slong nb_bad;
    acb_mat_t Ntau;
    acb_t scal;
    slong g = acb_mat_nrows(tau);
    slong n = 1<<g;
    slong i;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    
    acb_mat_init(Ntau, g, g);
    acb_init(scal);

    fmpz_mat_set(acb_theta_agm_ctx_matrix(ctx, k), N);
    acb_siegel_transform(Ntau, N, tau, prec);

    flint_printf("(ctx_set_const)\n");
    fmpz_mat_print_pretty(N); flint_printf("\n");
    acb_mat_printd(Ntau, 10); flint_printf("\n");
    
    nb_bad = FLINT_MAX(1, acb_theta_agm_nb_bad_steps(Ntau, prec));
    acb_theta_agm_ctx_reset_steps(ctx, k, nb_bad);

    /* Set roots to low precision, correctly rescaled */
    for (i = 0; i < nb_bad; i++)
    {
        acb_theta_naive_const(acb_theta_agm_ctx_roots(ctx, k) + i*n,
                Ntau, lowprec);
        acb_mat_scalar_mul_2exp_si(Ntau, Ntau, 1);
    }
    acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], lowprec);
    _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
            acb_theta_agm_ctx_roots(ctx, k), n * nb_bad, scal, lowprec);
  
    acb_mat_clear(Ntau);
    acb_clear(scal);
}

void
acb_theta_agm_ctx_set_const(acb_theta_agm_ctx_t ctx, const acb_mat_t tau,
        slong prec)
{  
    acb_mat_t half;
    acb_ptr th;
    arf_t rad;
    slong lowprec = ACB_THETA_AGM_LOWPREC;  
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong g = acb_mat_nrows(tau);
    slong k;
    int res;
    slong try = -1;

    acb_mat_init(half, g, g);
    th = _acb_vec_init(1<<g);

    acb_mat_scalar_mul_2exp_si(half, tau, -1);
    acb_theta_naive_const_proj(th, half, prec);

    while (try < ACB_THETA_AGM_NB_MATRIX_SETUPS)
    {
        try++;
        acb_theta_agm_ctx_candidates(acb_theta_agm_ctx_matrix(ctx, 0), try, g);
        arf_pos_inf(acb_theta_agm_ctx_rho(ctx));
        arf_zero(acb_theta_agm_ctx_max(ctx));
      
        for (k = 0; k < n; k++)
	{
            set_matrix(ctx, k, tau, acb_theta_agm_ctx_matrix(ctx, k), prec);
            acb_theta_agm_ctx_update_bounds(ctx, k, lowprec);
      
            acb_theta_agm_radius(rad, acb_theta_agm_ctx_mi(ctx, k),
                    acb_theta_agm_ctx_M0(ctx, k),
                    acb_theta_agm_ctx_minf(ctx, k),
                    acb_theta_agm_ctx_nb_bad_steps(ctx, k), prec);
            acb_theta_dupl_transform_radius_const(rad, rad, th,
                    acb_theta_agm_ctx_matrix(ctx, k), lowprec);
            arf_min(acb_theta_agm_ctx_rho(ctx),
                    acb_theta_agm_ctx_rho(ctx), rad);
            
            arf_div(rad, acb_theta_agm_ctx_M0(ctx, k),
                    acb_theta_agm_ctx_minf(ctx, k), lowprec, ARF_RND_CEIL);
            arf_mul_2exp_si(rad, rad, 1);
            arf_max(acb_theta_agm_ctx_max(ctx),
                    acb_theta_agm_ctx_max(ctx), rad);
	}

        if (!arf_is_finite(acb_theta_agm_ctx_max(ctx))) continue;
        if (arf_cmp_si(acb_theta_agm_ctx_rho(ctx), 0) <= 0) continue;

        res = acb_theta_agm_ctx_set_inv_der(ctx, th, prec);
        if (res) break;
    }

    acb_mat_clear(half);
    _acb_vec_clear(th, 1<<g);
}
