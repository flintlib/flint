
#include "acb_theta.h"


/* Compute the target of Newton scheme to high precision */

static void
acb_theta_newton_target(acb_ptr im, const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong dim = acb_theta_agm_ctx_dim(ctx);
    slong n = 1<<g;
    slong k;
    fmpz_t eps;
    acb_t zeta;
    
    fmpz_init(eps);
    acb_init(zeta);

    for (k = 0; k < n-1; k++)
    {
        acb_onei(zeta);  
        acb_pow_fmpz(zeta, zeta, acb_theta_agm_ctx_eps(ctx, k+1), prec);

        if (acb_theta_agm_ctx_is_ext(ctx))
        {
            acb_theta_transform_scal(&im[k], &im[k+n-1],
                    acb_theta_agm_ctx_z(ctx),
                    acb_theta_agm_ctx_tau(ctx),
                    acb_theta_agm_ctx_mat(ctx, k+1),
                    acb_theta_agm_ctx_k2(ctx, k+1), prec);
            acb_mul(&im[k], &im[k], zeta, prec);
            acb_mul(&im[k+n-1], &im[k+n-1], zeta, prec);
        }
        else
        {            
            acb_theta_transform_scal_const(&im[k], acb_theta_agm_ctx_tau(ctx),
                    acb_theta_agm_ctx_mat(ctx, k+1),
                    acb_theta_agm_ctx_k2(ctx, k+1), prec);
            acb_mul(&im[k], &im[k], zeta, prec);
        }
    }
    for (k = 0; k < dim; k++) acb_inv(&im[k], &im[k], prec);
    
    flint_printf("Target:\n");
    for (k = 0; k < dim; k++)
    {
        acb_printd(&im[k], 10); flint_printf("\n");
    }
    
    fmpz_clear(eps);
    acb_clear(zeta);
}

/* Start Newton scheme. Output: im, start are exact. Return the absolute
   precision of start. */

static slong
acb_theta_newton_start(acb_ptr start, acb_ptr im, arf_t err,
        const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong dim = acb_theta_agm_ctx_dim(ctx);
    slong n = 1<<g;
    slong log_M, log_B2, log_B3;
    arf_t e;
    fmpz_t exp;
    acb_mat_t half;
    slong k;

    if (acb_theta_agm_ctx_is_ext(ctx)) n *= 2;
    arf_init(e);
    acb_mat_init(half, g, g);
    log_M = acb_theta_agm_ctx_log_M(ctx);
    log_B2 = acb_theta_agm_ctx_log_B2(ctx);
    log_B3 = acb_theta_agm_ctx_log_B3(ctx);
    
    acb_mat_scalar_mul_2exp_si(half, acb_theta_agm_ctx_tau(ctx), -1);
  
    /* Get image; add some error; get midpoint and error */
    acb_theta_newton_target(im, ctx, prec);
    arf_one(e);
    arf_mul_2exp_si(e, e, -prec-1-log_B3);
    for (k = 0; k < dim; k++) acb_add_error_arf(&im[k], e);
  
    arf_zero(err);
    for (k = 0; k < dim; k++)
    {
        arf_set_mag(e, arb_radref(acb_realref(&im[k])));
        arf_max(err, err, e);
        arf_set_mag(e, arb_radref(acb_imagref(&im[k])));
        arf_max(err, err, e);
        acb_get_mid(&im[k], &im[k]);
    }
    arf_mul_2exp_si(err, err, 1);
    arf_frexp(e, exp, err);
    prec = -fmpz_get_si(exp);

    flint_printf("newton_start: chosen target prec %wd\n", prec);
  
    /* im is now exact, and known to precision prec. Pick starting precision */
    while ((prec > ACB_THETA_AGM_BASEPREC - log_M - ACB_THETA_AGM_GUARD)
            && (prec > 2*(log_B2 + log_B3 + 2)))
    {
        prec = (prec + log_B2 + log_B3 + 3)/2;
    }
    flint_printf("newton_start: starting prec %wd\n", prec);

    /* Set start using naive algorithm; control error bound; get midpoints */
    _acb_vec_set(start, acb_theta_agm_ctx_th(ctx), n);    
    for (k = 0; k < n; k++)
    {
        if (mag_cmp_2exp_si(arb_radref(acb_realref(&start[k])), -prec-1) > 0
                || mag_cmp_2exp_si(arb_radref(acb_imagref(&start[k])), -prec-1) > 0)
	{
            flint_printf("acb_theta_newton_start: Error (insufficient precision)\n");
            fflush(stdout);
            flint_abort();
	}
        acb_get_mid(&start[k], &start[k]);
    }
  
    arf_clear(e);
    acb_mat_clear(half);
    return prec;
}

/* Newton step. Input: current is exact, and an approximation of desired output
   to absolute precision 2^-prec. im is exact. Output: current is exact, and
   an approximation of desired output to a superior precision (return value) */

static slong
acb_theta_newton_step(acb_ptr next, acb_srcptr current, acb_srcptr im,
        const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong dim = acb_theta_agm_ctx_dim(ctx);
    slong n = 1<<g;
    slong log_th, log_M, log_B1, log_B2, log_B3;
    slong log_eta, nextprec;
    arb_t eta;
    acb_mat_t fd;
    acb_ptr f;
    acb_mat_t h;
    int res;
    slong k;

    arb_init(eta);
    acb_mat_init(fd, dim, dim);
    f = _acb_vec_init(dim);
    acb_mat_init(h, dim, 1);
    
    log_M = acb_theta_agm_ctx_log_M(ctx);
    log_B1 = acb_theta_agm_ctx_log_B1(ctx);
    log_B2 = acb_theta_agm_ctx_log_B2(ctx);
    log_B3 = acb_theta_agm_ctx_log_B3(ctx);
    log_th = acb_theta_agm_ctx_log_th(ctx);

    /* Set nextprec, eta */
    nextprec = 2*prec + 2*n_clog(dim, 2) + 2*log_B1 + 2*log_B3 + 9 + log_M
        + ACB_THETA_AGM_GUARD;
    log_eta = -(prec + log_B1 + log_B3 + n_clog(dim, 2) + 2);
    arb_one(eta);
    arb_mul_2exp_si(eta, eta, log_eta);

    flint_printf("(newton_step) log Bi: %wd %wd %wd\n", log_B1, log_B2, log_B3);
    flint_printf("log_rho, log_th, log_M: %wd %wd %wd\n",
            acb_theta_agm_ctx_log_rho(ctx), log_th, log_M);

    /* Compute correction */
    acb_theta_newton_fd(f, fd, current, eta, ctx, nextprec);
    
    flint_printf("Current image:\n");
    for (k = 0; k < dim; k++)
    {
        acb_printd(&f[k], 10); flint_printf("\n");
    }
    
    res = acb_mat_inv(fd, fd, nextprec);
    if (!res)
    {
        flint_printf("acb_theta_newton_step: Error (impossible inversion)\n");
        fflush(stdout);
        flint_abort();
    }
    _acb_vec_sub(f, im, f, dim, nextprec);
    
    flint_printf("Current error:\n");
    for (k = 0; k < dim; k++)
    {
        acb_printd(&f[k], 10); flint_printf("\n");
    }

    for (k = 0; k < dim; k++)
    {
        acb_set(acb_mat_entry(h, k, 0), &f[k]);
    }
    acb_mat_mul(h, fd, h, nextprec);
  
    /* Check that h does not have too much additional error */
    nextprec = 2*prec - log_B2 - log_B3 - 2;
    for (k = 0; k < dim; k++)
    {
        if (mag_cmp_2exp_si(
                        arb_radref(acb_realref(acb_mat_entry(h, k, 0))),
                        -nextprec-1) > 0
                || mag_cmp_2exp_si(
                        arb_radref(acb_imagref(acb_mat_entry(h, k, 0))),
                        -nextprec-1) > 0)
	{
            flint_printf("acb_theta_newton_step: Error (imprecise correction)\n");
            flint_printf("Needed prec %wd\n", nextprec+1);
            fflush(stdout);
            flint_abort();
	}
    }

    /* Set result */
    for (k = 0; k < n; k++)
    {
        acb_set(&next[k], &current[k]);
        if (k > 0)
        {
            acb_add(&next[k], &next[k], acb_mat_entry(h, k-1, 0),
                nextprec + log_th + ACB_THETA_AGM_GUARD);
        }
        acb_get_mid(&next[k], &next[k]);
        
        if (acb_theta_agm_ctx_is_ext(ctx))
        {
            if (k > 0)
            {                
                acb_add(&next[k+n], &next[k+n], acb_mat_entry(h, k+n-2, 0),
                        nextprec + log_th + ACB_THETA_AGM_GUARD);
            }
            acb_get_mid(&next[k+n], &next[k+n]);
        }
    }
    flint_printf("Precision increase from %wd to %wd\n", prec, nextprec);
  
    arb_clear(eta);
    acb_mat_clear(fd);
    _acb_vec_clear(f, dim);
    acb_mat_clear(h);  
    return nextprec;
}


void
acb_theta_newton_run(acb_ptr r, const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = 1<<g;
    slong dim = acb_theta_agm_ctx_dim(ctx);
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);
    slong log_B3;
    acb_ptr im;
    arf_t err;
    fmpz_t exp;
    slong current_prec;
    slong k;

    im = _acb_vec_init(dim);
    arf_init(err);
    fmpz_init(exp);

    log_B3 = acb_theta_agm_ctx_log_B3(ctx);
    current_prec = acb_theta_newton_start(r, im, err, ctx, prec);
    arf_frexp(err, exp, err);
    prec = -fmpz_get_si(exp);
    
    while (current_prec < prec)
    {
        current_prec = acb_theta_newton_step(r, r, im, ctx, current_prec);
    }
    /* Add error: coming from prec, and coming from err */
    arf_one(err);
    arf_mul_2exp_si(err, err, -current_prec + log_B3 + 1);

    flint_printf("(newton_run) Before error\n");
    for (k = 0; k < n; k++)
    {
        acb_printd(&r[k], 10); flint_printf("\n");        
    }
    if (is_ext)
    {
        for (k = 0; k < n; k++)
        {
            acb_printd(&r[k+n], 10); flint_printf("\n");        
        }
    }
        

    flint_printf("(newton_run) Additional error\n");
    arf_printd(err, 10); flint_printf("\n");
    
    for (k = 1; k < n; k++)
    {
        acb_add_error_arf(&r[k], err);
        if (is_ext) acb_add_error_arf(&r[k+n], err);
    }
    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);
    for (k = 1; k < n; k++)
    {
        acb_add_error_arf(&r[k], err);
        if (is_ext) acb_add_error_arf(&r[k+n], err);
    }

    _acb_vec_clear(im, dim);
    arf_clear(err);
    fmpz_clear(exp);
}
