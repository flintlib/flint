
#include "acb_theta.h"

/* Collect data for a given symplectic matrix */

static void
set_matrix(acb_theta_agm_ctx_t ctx, slong k, acb_srcptr z,
        const acb_mat_t tau, const fmpz_mat_t N, slong prec)
{
    slong nb_bad;
    acb_mat_t Ntau;
    acb_ptr Nz;
    acb_t scal;
    slong g = acb_mat_nrows(tau);
    slong n = 1<<(g+1);
    slong i;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
  
    acb_mat_init(Ntau, g, g);
    Nz = _acb_vec_init(g);
    acb_init(scal);

    fmpz_mat_set(acb_theta_agm_ctx_matrix(ctx, k), N);
    acb_siegel_transform_ext(Nz, Ntau, N, z, tau, prec);

    flint_printf("(ctx_set)\n");
    fmpz_mat_print_pretty(N); flint_printf("\n");
    acb_mat_printd(Ntau, 10); flint_printf("\n");
    for (i = 0; i < g; i++)
    {
        acb_printd(&Nz[i], 10); flint_printf("\n");
    }
    
    nb_bad = FLINT_MAX(1, acb_theta_agm_ext_nb_bad_steps(Nz, Ntau, prec));
    acb_theta_agm_ctx_reset_steps(ctx, k, nb_bad);

    /* Set roots to low precision, correctly rescaled */
    for (i = 0; i < nb_bad; i++)
    {
        acb_theta_naive(acb_theta_agm_ctx_roots(ctx, k) + i*n,
                Nz, 1, Ntau, lowprec);
        acb_mat_scalar_mul_2exp_si(Ntau, Ntau, 1);
    }
    acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], lowprec);
    _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
            acb_theta_agm_ctx_roots(ctx, k), n * nb_bad, scal, lowprec);

    acb_mat_clear(Ntau);
    _acb_vec_clear(Nz, g);
    acb_clear(scal);
}

void
acb_theta_agm_ctx_set(acb_theta_agm_ctx_t ctx, acb_srcptr z,
        const acb_mat_t tau, slong prec)
{  
    acb_mat_t half;
    acb_ptr th;
    arf_t rad;
    arb_t max;
    slong lowprec = ACB_THETA_AGM_LOWPREC;  
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong g = acb_mat_nrows(tau);
    slong k;
    int res;
    slong try = -1;

    acb_mat_init(half, g, g);
    th = _acb_vec_init(1<<(g+1));
    arb_init(max);
    arf_init(rad);

    acb_mat_scalar_mul_2exp_si(half, tau, -1);
    acb_theta_naive_proj(th, z, 1, half, prec);

    while (try < ACB_THETA_AGM_NB_MATRIX_SETUPS)
    {
        try++;        
        acb_theta_agm_ctx_candidates(acb_theta_agm_ctx_matrix(ctx, 0), try, g);
        arf_pos_inf(acb_theta_agm_ctx_rho(ctx));
        arf_zero(acb_theta_agm_ctx_max(ctx));
      
        for (k = 0; k < n; k++)
	{
            set_matrix(ctx, k, z, tau, acb_theta_agm_ctx_matrix(ctx, k), prec);
            acb_theta_agm_ctx_update_bounds(ctx, k, lowprec);
         
            acb_theta_agm_radius(rad, acb_theta_agm_ctx_mi(ctx, k),
                    acb_theta_agm_ctx_M0(ctx, k),
                    acb_theta_agm_ctx_minf(ctx, k),
                    acb_theta_agm_ctx_nb_bad_steps(ctx, k), lowprec);
            acb_theta_dupl_transform_radius(rad, rad, th,
                    acb_theta_agm_ctx_matrix(ctx, k), lowprec);
            arf_min(acb_theta_agm_ctx_rho(ctx),
                    acb_theta_agm_ctx_rho(ctx), rad);

            arb_set_arf(max, acb_theta_agm_ctx_M0(ctx, k));
            arb_add_arf(max, max, rad, lowprec);
            arb_div_arf(max, max, acb_theta_agm_ctx_minf(ctx, k), lowprec);
            arb_mul_2exp_si(max, max, 3);
            arb_log(max, max, lowprec);
            arb_sqr(max, max, lowprec);
            arb_mul_si(max, max, 20, lowprec);
            arb_exp(max, max, lowprec);
            arb_get_ubound_arf(rad, max, lowprec);
            arf_max(acb_theta_agm_ctx_max(ctx),
                    acb_theta_agm_ctx_max(ctx), rad);            
	}

        if (!arf_is_finite(acb_theta_agm_ctx_max(ctx))) continue;
        if (arf_cmp_si(acb_theta_agm_ctx_rho(ctx), 0) <= 0) continue;

        res = acb_theta_agm_ctx_set_inv_der(ctx, th, prec);
        if (res) break;
    }
    acb_mat_clear(half);
    _acb_vec_clear(th, 1<<(g+1));
    arb_clear(max);
    arf_clear(rad);
}
