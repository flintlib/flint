
#include "acb_theta.h"

void acb_theta_newton_run(acb_ptr r, const acb_mat_t tau, const acb_theta_agm_ctx_t ctx,
			  slong prec)
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
