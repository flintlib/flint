
#include "acb_theta.h"

/* Output: im, start are exact. Return the absolute precision of start. */

static void acb_theta_newton_target(acb_ptr im, const acb_mat_t tau,
				    const acb_theta_agm_ctx_t ctx, slong prec);
{
  slong g = acb_theta_agm_ctx_g(ctx);
  slong n = acb_theta_agm_ctx_nb(ctx);
  slong k;
  acb_mat_t w;
  fmpz_t epsilon;
  acb_t zeta8, mu;

  acb_mat_init(w, g, g);
  fmpz_init(epsilon);
  acb_init(zeta8);
  acb_init(mu);

  acb_one(zeta);
  acb_mul_2exp_si(zeta, zeta, -2);
  acb_exp_pi_i(zeta, zeta, prec);
  
  for (k = 0; k < n; k++)
    {
      acb_theta_transform_image_char(epsilon, 0, acb_theta_agm_ctx_matrix(ctx, k));
      acb_pow_si(mu, zeta, fmpz_get_si(epsilon), prec);
      acb_siegel_cocycle(w, acb_theta_agm_ctx_matrix(ctx, k), tau, prec);
      acb_mat_det(&im[k], w, prec);
      acb_mul(&im[k], &im[k], mu, prec);
    }

  acb_mat_clear(w);
  fmpz_clear(epsilon);
  acb_clear(zeta8);
  acb_clear(mu);
}

slong acb_theta_newton_start(acb_ptr start, acb_ptr im, arf_t err, const acb_mat_t tau,
			     const acb_theta_agm_ctx_t ctx, slong prec)
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
  acb_theta_naive_const_proj(start, half, prec + log_max + ACB_THETA_AGM_GUARD);
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
}
