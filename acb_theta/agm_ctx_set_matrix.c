
#include "acb_theta.h"

void
acb_theta_agm_ctx_set_matrix(acb_theta_agm_ctx_t ctx, slong k,
        const acb_mat_t tau, const fmpz_mat_t N, slong prec)
{
    slong nb_bad;
    acb_mat_t z;
    acb_t scal;
    arb_t abs;
    arb_t m;
    slong g = acb_mat_nrows(tau);
    slong n = 1<<g;
    slong i, j;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
  
    acb_mat_init(z, g, g);
    arb_init(abs);
    arb_init(m);
    acb_init(scal);

    fmpz_mat_set(acb_theta_agm_ctx_matrix(ctx, k), N);

    /* Compute number of bad steps */
    acb_siegel_transform(z, N, tau, prec);
    nb_bad = FLINT_MIN(1, acb_theta_agm_nb_bad_steps(z, prec));
    acb_theta_agm_ctx_reset_steps(ctx, k, nb_bad);

    /* Set roots to low precision, correctly rescaled */
    for (i = 0; i < nb_bad; i++)
    {
        acb_mat_printd(z, 10); flint_printf("\n");
        acb_theta_naive_const(acb_theta_agm_ctx_roots(ctx, k) + i*n,
                z, lowprec);
        acb_mat_scalar_mul_2exp_si(z, z, 1);
    }
    acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], lowprec);
    _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
            acb_theta_agm_ctx_roots(ctx, k), n * nb_bad, scal, lowprec);

    /* Compute M0 */
    arb_zero(m);
    for (i = 0; i < n*nb_bad; i++)
    {
        acb_abs(abs, &acb_theta_agm_ctx_roots(ctx, k)[i], lowprec);
        arb_max(m, m, abs, lowprec);
    }
    arb_sqr(m, m, lowprec);
    arb_get_ubound_arf(acb_theta_agm_ctx_M0(ctx, k), m, lowprec);

    /* Compute minf */
    acb_abs(m, &acb_theta_agm_ctx_roots(ctx, k)[n*(nb_bad-1)], lowprec);
    arb_div_si(m, m, 20, lowprec); /* Cf. agm_nb_bad_steps */
    arb_get_lbound_arf(acb_theta_agm_ctx_minf(ctx, k), m, lowprec);

    /* Compute mi */
    for (i = 0; i < nb_bad; i++)
    {
        arb_pos_inf(m);
        for (j = 0; j < n; j++)
	{
            acb_abs(abs, &acb_theta_agm_ctx_roots(ctx,k)[i*n + j], lowprec);
            arb_min(m, m, abs, lowprec);
	}
        arb_sqr(m, m, lowprec);
        arb_get_lbound_arf(&acb_theta_agm_ctx_mi(ctx,k)[i], m, lowprec);
    }
  
    acb_mat_clear(z);
    acb_clear(scal);
    arb_clear(abs);
    arb_clear(m);
}
