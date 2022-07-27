
#include "acb_theta.h"

/* We implement rigorous Newton iterations.

   Input: current is exact, and an approximation of desired output to
   absolute precision 2^-prec. im is exact.

   Output: current is exact, and an approximation of desired output to
   a superior precision (as in paper) */

slong acb_theta_newton_step(acb_ptr next, acb_srcptr current, acb_srcptr im,
			    const acb_theta_agm_ctx_t, slong prec)
{
  slong g = acb_theta_agm_ctx_g(ctx);
  slong n = acb_theta_agm_ctx_nb(ctx); /* dimension is n-1 */
  slong log_max, log_rho, log_B1, log_B2, log_B3;
  slong log_eta, nextprec, nprime;
  arb_t eta;
  acb_mat_t fd;
  acb_ptr f;
  acb_mat_t h;
  int res;

  arb_init(eta);
  acb_mat_init(fd, n-1, n-1);
  f = _acb_vec_init(n-1);
  acb_mat_init(h, n-1, 1);

  /* Set logs */
  acb_theta_newton_logs(&log_max, &log_rho, &log_B1, &log_B2, &log_B3, ctx);

  /* Set nextprec, eta */
  nextprec = 2*prec + 2*n_clog(n,2) + 2*log_B1 + 2*log_B3 + 9 + log_max + ACB_THETA_AGM_GUARD;
  log_eta = -(prec + log_B1 + log_B3 + n_clog(n, 2) + 2);
  arb_one(eta);
  arb_mul_2exp_si(eta, eta, log_eta);

  /* Compute correction */
  acb_theta_newton_fd(r, fd, current, eta, ctx, nextprec);
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
      if (mag_cmp_2exp_si(arb_magref(acb_realref(acb_mat_entry(h, k, 0))), -nprime-1) > 0
	  || mag_cmp_2exp_si(arb_magref(acb_imagref(acb_mat_entry(h, k, 0))), -nprime-1) > 0)
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
