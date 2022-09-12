
#include "acb_theta.h"

/* Compute radius of disk around (theta_i^2/theta_0^2(N.tau)) where
   AGM function is surely well-defined */

static void
agm_radius(arf_t rad, const arf_struct* mi, const arf_t M0,
        const arf_t minf, slong nb, slong prec)
{
    arf_t prod, term, res;
    slong j;

    arf_init(prod);
    arf_init(term);
    arf_init(res);
  
    arf_one(prod);
    arf_mul_2exp_si(res, &mi[0], -1);	  
    for (j = 0; j < nb; j++)
    {
        arf_mul_2exp_si(term, M0, 1);
        arf_add(term, term, &mi[j], prec, ARF_RND_CEIL);
        arf_div(term, &mi[j], term, prec, ARF_RND_FLOOR);
        arf_sqrt(term, term, prec, ARF_RND_FLOOR);
        arf_mul(prod, prod, term, prec, ARF_RND_FLOOR);
      
        if (j == nb - 1) arf_mul(term, minf, prod, prec, ARF_RND_FLOOR);
        else arf_mul(term, &mi[j+1], prod, prec, ARF_RND_FLOOR);
        arf_mul_2exp_si(term, term, -1);
        arf_min(res, res, term);
    }

    arf_set(rad, res);
    arf_clear(prod);
    arf_clear(term);
    arf_clear(res);
}

/* Given result r of agm_radius, compute radius of disk around
   (theta_i/theta_0(tau/2)) where dupl+transform+agm is surely
   well-defined */

static void
propagate_rho(arf_t rho, const arf_t r, acb_srcptr th_proj,
        const fmpz_mat_t N, slong prec)
{
    ulong ab_0, ab;
    fmpz_t epsilon;
    acb_ptr th_dupl;
    arb_t abs_0, abs;
    arf_t bound, max, res;
    slong g = fmpz_mat_nrows(N)/2;
    slong k;

    fmpz_init(epsilon);
    th_dupl = _acb_vec_init(1<<(2*g));
    arb_init(abs_0);
    arb_init(abs);
    arf_init(bound);
    arf_init(max);
    arf_init(res);

    acb_theta_duplication_all(th_dupl, th_proj, g, prec);
  
    /* theta_i^2/theta_0^2 is obtained as quotient of two theta
       transformations, multiplied by some root of unity */
    ab_0 = acb_theta_transform_image_char(epsilon, 0, N);
    acb_abs(abs_0, &th_dupl[ab_0], prec);
    arf_pos_inf(res);
  
    for (k = 1; k < (1<<g); k++)
    {
        ab = acb_theta_transform_image_char(epsilon, 0, N);
        acb_abs(abs, &th_dupl[ab], prec);
        arb_add(abs, abs, abs_0, prec);
        arb_div(abs, abs_0, abs, prec);
        arb_mul(abs, abs, abs_0, prec);
        arb_min(abs, abs, abs_0, prec);
        arb_mul_2exp_si(abs, abs, -1);
        arb_get_lbound_arf(bound, abs, prec);
        arf_min(res, res, bound);
    }

    arf_zero(max);
    /* Now we have a suitable radius for the duplicated values */
    for (k = 1; k < (1<<g); k++)
    {
        acb_abs(abs, &th_proj[k], prec);
        arb_get_ubound_arf(bound, abs, prec);
        arf_max(max, max, bound);
    }
    arf_div(res, res, max, prec, ARF_RND_FLOOR);
    arf_div_si(res, res, 3, prec, ARF_RND_FLOOR);
    arf_min(res, res, max);

    fmpz_clear(epsilon);
    _acb_vec_clear(th_dupl, 1<<(2*g));
    arb_clear(abs_0);
    arb_clear(abs);
    arf_clear(bound);
    arf_clear(max);
    arf_clear(res);  
}


void
acb_theta_agm_ctx_set_all(acb_theta_agm_ctx_t ctx, const acb_mat_t tau,
        slong prec)
{  
    acb_mat_t half;
    acb_ptr th;
    arf_t rad;
    arf_t B2;
    fmpz_t e;
    slong exp;
    arb_t eta;
    acb_ptr r;
    acb_mat_t fd, fdinv;
    arb_t norm, bound, test;  
    slong lowprec = ACB_THETA_AGM_LOWPREC;  
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong g = acb_mat_nrows(tau);
    slong k;
    int stop = 0;
    int res;
    int try = 0;

    acb_mat_init(half, g, g);
    th = _acb_vec_init(1<<g);
    arf_init(rad);
    arf_init(B2);
    fmpz_init(e);
    arb_init(eta);
    r = _acb_vec_init(n-1);
    acb_mat_init(fd, n-1, n-1);
    acb_mat_init(fdinv, n-1, n-1);
    arb_init(norm);
    arb_init(bound);
    arb_init(test);

    acb_mat_scalar_mul_2exp_si(half, tau, -1);
    acb_theta_naive_const_proj(th, half, prec);

    while (!stop && (try < ACB_THETA_AGM_NB_MATRIX_SETUPS))
    {
        try++;
        acb_theta_agm_ctx_matrices(acb_theta_agm_ctx_matrix(ctx, 0), try, g);
        arf_pos_inf(acb_theta_agm_ctx_rho(ctx));
        arf_zero(acb_theta_agm_ctx_max(ctx));
      
        for (k = 0; k < n; k++)
	{
            acb_theta_agm_ctx_set_matrix(ctx, k, tau,
                    acb_theta_agm_ctx_matrix(ctx, k), prec);
	}
      
        for (k = 0; k < n; k++)
	{
            agm_radius(rad, acb_theta_agm_ctx_mi(ctx, k),
                    acb_theta_agm_ctx_M0(ctx, k),
                    acb_theta_agm_ctx_minf(ctx, k),
                    acb_theta_agm_ctx_nb_bad_steps(ctx, k), prec);
	  
            /* Propagate radius according to quotients & duplication */
            propagate_rho(rad, rad, th, acb_theta_agm_ctx_matrix(ctx, k),
                    lowprec);
            arf_min(acb_theta_agm_ctx_rho(ctx),
                    acb_theta_agm_ctx_rho(ctx), rad);
	  
            /* Update maximum value of Borchardt quotients */
            arf_div(rad, acb_theta_agm_ctx_M0(ctx, k),
                    acb_theta_agm_ctx_minf(ctx, k), lowprec, ARF_RND_CEIL);
            arf_mul_2exp_si(rad, rad, 1);
            arf_max(acb_theta_agm_ctx_max(ctx),
                    acb_theta_agm_ctx_max(ctx), rad);
	}

        if (!arf_is_finite(acb_theta_agm_ctx_max(ctx))) continue;
        if (arf_cmp_si(acb_theta_agm_ctx_rho(ctx), 0) <= 0) continue;

        /* Evaluate finite difference */
        acb_theta_cauchy(B2, acb_theta_agm_ctx_rho(ctx),
                acb_theta_agm_ctx_max(ctx), 2, n-1, lowprec);
        arf_frexp(B2, e, B2);
        exp = fmpz_get_si(e);
        arb_one(eta);
        arb_mul_2exp_si(eta, eta, FLINT_MIN(-exp-n_clog(n,2), -prec/2));
        acb_theta_newton_fd(r, fd, th, eta, ctx, prec);
        res = acb_mat_inv(fdinv, fd, prec);
        if (!res) continue;
      
        /* Is ||FD^-1||*n*B2*eta less than 1? */
        acb_mat_ninf(norm, fdinv, lowprec);
        arb_mul_arf(bound, norm, B2, lowprec);
        arb_mul_si(bound, bound, n, lowprec);
        arb_mul(bound, bound, eta, lowprec);
        arb_sub_si(test, bound, 1, lowprec);
        if (!arb_is_negative(test)) continue;

        /* Get inv_der */
        arb_mul(bound, bound, norm, lowprec);
        arb_add(bound, bound, norm, lowprec);
        arb_get_ubound_arf(acb_theta_agm_ctx_inv_der(ctx), bound, lowprec);
        stop = 1;
    }

    /* Clear */
    acb_mat_clear(half);
    _acb_vec_clear(th, 1<<g);
    arf_clear(rad);
    arf_clear(B2);
    fmpz_clear(e);
    arb_clear(eta);
    _acb_vec_clear(r, n-1);
    acb_mat_clear(fd);
    acb_mat_clear(fdinv);
    arb_clear(norm);
    arb_clear(bound);
    arb_clear(test);
}
