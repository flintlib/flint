
#include "acb_theta.h"

/* Get binary logs of data stored in ctx */

static void
acb_theta_newton_logs(slong* log_max, slong* log_rho, slong* log_B1,
        slong* log_B2, slong* log_B3, const acb_theta_agm_ctx_t ctx)
{
    arf_t c;
    fmpz_t e;
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong lowprec = ACB_THETA_AGM_LOWPREC;

    arf_init(c);
    fmpz_init(e);
  
    arf_frexp(c, e, acb_theta_agm_ctx_max(ctx));
    *log_max = fmpz_get_si(e);
    arf_frexp(c, e, acb_theta_agm_ctx_rho(ctx));
    *log_rho = fmpz_get_si(e) - 1;  
    arf_mul_si(c, acb_theta_agm_ctx_max(ctx), 2*n, lowprec, ARF_RND_CEIL);
    arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
    arf_frexp(c, e, c);
    *log_B1 = fmpz_get_si(e);
    arf_mul_si(c, acb_theta_agm_ctx_max(ctx), 2*n*(n+1), lowprec,
            ARF_RND_CEIL);
    arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
    arf_div(c, c, acb_theta_agm_ctx_rho(ctx), lowprec, ARF_RND_CEIL);
    arf_frexp(c, e, c);
    *log_B2 = fmpz_get_si(e);  
    arf_frexp(c, e, acb_theta_agm_ctx_inv_der(ctx));
    *log_B3 = fmpz_get_si(e);

    arf_clear(c);
    fmpz_clear(e);
}

/* Compute the target of Newton scheme to low precision */

static void
acb_theta_newton_target(acb_ptr im, const acb_mat_t tau,
        const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong k;
    acb_mat_t w;
    fmpz_t epsilon;
    acb_t zeta, mu;

    acb_mat_init(w, g, g);
    fmpz_init(epsilon);
    acb_init(zeta);
    acb_init(mu);

    acb_one(zeta);
    acb_mul_2exp_si(zeta, zeta, -2);
    acb_exp_pi_i(zeta, zeta, prec);
  
    for (k = 0; k < n; k++)
    {
        acb_theta_transform_image_char(epsilon, 0,
                acb_theta_agm_ctx_matrix(ctx, k));
        acb_pow_si(mu, zeta, fmpz_get_si(epsilon), prec);
        acb_siegel_cocycle(w, acb_theta_agm_ctx_matrix(ctx, k), tau, prec);
        acb_mat_det(&im[k], w, prec);
        acb_mul(&im[k], &im[k], mu, prec);
    }

    acb_mat_clear(w);
    fmpz_clear(epsilon);
    acb_clear(zeta);
    acb_clear(mu);
}

/* Start Newton scheme. Output: im, start are exact. Return the absolute
   precision of start. */

static slong
acb_theta_newton_start(acb_ptr start, acb_ptr im, arf_t err,
        const acb_mat_t tau, const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong log_max, log_rho, log_B1, log_B2, log_B3;
    arf_t e;
    fmpz_t exp;
    acb_mat_t half;
    slong k;

    arf_init(e);
    acb_mat_init(half, g, g);
  
    acb_theta_newton_logs(&log_max, &log_rho, &log_B1, &log_B2, &log_B3, ctx);
    acb_mat_scalar_mul_2exp_si(half, tau, -1);
  
    /* Get image; add some error; get midpoint and error */
    acb_theta_newton_target(im, tau, ctx, prec);
    arf_one(e);
    arf_mul_2exp_si(e, e, -prec-1-log_B3);
    for (k = 0; k < n-1; k++) acb_add_error_arf(&im[k], e);
  
    arf_zero(err);
    for (k = 0; k < n-1; k++)
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
  
    /* im is now exact, and known to precision prec. Pick starting precision */
    while ((prec > ACB_THETA_AGM_BASEPREC)
            && (prec > 2*(log_B2 + log_B3 + 2)))
    {
        prec = (prec + log_B2 + log_B3 + 3)/2;
    }

    /* Set start using naive algorithm; control error bound; get midpoints */
    acb_theta_naive_const_proj(start, half,
            prec + log_max + ACB_THETA_AGM_GUARD);
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
    slong n = acb_theta_agm_ctx_nb(ctx); /* dimension is n-1 */
    slong log_max, log_rho, log_B1, log_B2, log_B3;
    slong log_eta, nextprec, nprime;
    arb_t eta;
    acb_mat_t fd;
    acb_ptr f;
    acb_mat_t h;
    int res;
    slong k;

    arb_init(eta);
    acb_mat_init(fd, n-1, n-1);
    f = _acb_vec_init(n-1);
    acb_mat_init(h, n-1, 1);

    /* Set logs */
    acb_theta_newton_logs(&log_max, &log_rho, &log_B1, &log_B2, &log_B3, ctx);

    /* Set nextprec, eta */
    nextprec = 2*prec + 2*n_clog(n,2) + 2*log_B1 + 2*log_B3 + 9 + log_max
        + ACB_THETA_AGM_GUARD;
    log_eta = -(prec + log_B1 + log_B3 + n_clog(n, 2) + 2);
    arb_one(eta);
    arb_mul_2exp_si(eta, eta, log_eta);

    /* Compute correction */
    acb_theta_newton_fd(f, fd, current, eta, ctx, nextprec);
    res = acb_mat_inv(fd, fd, nextprec);
    if (!res)
    {
        flint_printf("acb_theta_newton_step: Error (impossible inversion)\n");
        fflush(stdout);
        flint_abort();
    }
    _acb_vec_sub(f, im, f, n-1, prec);

    for (k = 0; k < n-1; k++)
    {
        acb_set(acb_mat_entry(h, k, 0), &f[k]);
    }
    acb_mat_mul(h, fd, h, nextprec);
  
    /* Check that h does not have too much additional error */
    nprime = 2*n - log_B2 - log_B3 - 2;
    for (k = 0; k < n-1; k++)
    {
        if (mag_cmp_2exp_si(arb_radref(acb_realref(acb_mat_entry(h, k, 0))), -nprime-1) > 0
                || mag_cmp_2exp_si(arb_radref(acb_imagref(acb_mat_entry(h, k, 0))), -nprime-1) > 0)
	{
            flint_printf("acb_theta_newton_step: Error (imprecise correction)\n");
            fflush(stdout);
            flint_abort();
	}
    }

    /* Set result */
    for (k = 0; k < n-1; k++)
    {
        acb_add(&next[k], &current[k], acb_mat_entry(h, k, 0),
                nprime + log_max + ACB_THETA_AGM_GUARD);
        acb_get_mid(&next[k], &next[k]);
    }
  
    arb_clear(eta);
    acb_mat_clear(fd);
    _acb_vec_clear(f, n-1);
    acb_mat_clear(h);  
    return nprime;
}


void
acb_theta_newton_run(acb_ptr r, const acb_mat_t tau,
        const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong n = acb_theta_agm_ctx_nb(ctx);
    slong log_max, log_rho, log_B1, log_B2, log_B3;
    acb_ptr im;
    arf_t err;
    fmpz_t exp;
    slong current_prec;
    slong k;

    im = _acb_vec_init(n-1);
    arf_init(err);
    fmpz_init(exp);

    acb_theta_newton_logs(&log_max, &log_rho, &log_B1, &log_B2, &log_B3, ctx);
    current_prec = acb_theta_newton_start(r, im, err, tau, ctx, prec);
    while (current_prec < prec)
    {
        current_prec = acb_theta_newton_step(r, r, im, ctx, prec);
    }
    /* Add error: coming from prec, and coming from err */
    arf_frexp(err, exp, err);
    arf_one(err);
    arf_mul_2exp_si(err, err, - fmpz_get_si(exp) + log_B3 + 1);
    for (k = 0; k < n-1; k++) acb_add_error_arf(&r[k], err);
    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);
    for (k = 0; k < n-1; k++) acb_add_error_arf(&r[k], err);

    _acb_vec_clear(im, n-1);
    arf_clear(err);
    fmpz_clear(exp);
}
