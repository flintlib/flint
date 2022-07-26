
#include "acb_theta.h"

static void agm_radius(arf_t rad, const arf_struct* mi, slong nb, slong prec)
{	  
  arb_one(prod);
  arb_mul_2exp_si(rad, &mi[0], -1);	  
  for (j = 0; j < nb_bad_steps[k]; j++)
    {
      arb_mul_2exp_si(term, M0, 1);
      arb_add(term, term, &mi[j], lowprec);
      arb_div(term, &mi[j], term, lowprec);
      arb_sqrt(term, term, lowprec);
      arb_mul(prod, prod, term, lowprec);
      
      if (j == nb_bad_steps[k] - 1) arb_mul(term, minf, prod, lowprec);
      else arb_mul(term, &mi[j+1], prod, lowprec);
      arb_mul_2exp_si(term, term, -1);
      arb_min(rad, rad, term, lowprec);
    }
}

static void propagate_rho(arb_t rho, const arb_t r, acb_srcptr th_half,
			  const fmpz_mat_t N, slong prec);


void acb_theta_newton_set_all(acb_theta_newton_t ctx, const acb_mat_t tau, slong prec)
{  
  acb_mat_t half;
  acb_ptr th;
  arf_t rad;
  arf_t B2;
  fmpz_t e;
  slong exp;
  arb_t eta;
  acb_mat_t fd, fdinv;
  arb_t norm, bound, test;  
  slong lowprec = ACB_THETA_AGM_LOWPREC;  
  slong n = acb_theta_newton_nb(ctx);
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
  acb_mat_init(fd, n-1, n-1);
  acb_mat_init(fdinv, n-1, n-1);
  arb_init(norm);
  arb_init(bound);
  arb_init(test);

  acb_mat_scalar_mul_2exp_si(half, tau, -1);
  acb_theta_naive_const(th, half, prec);

  while (!stop && (try < ACB_THETA_AGM_NB_MATRIX_SETUPS))
    {
      try++;
      acb_theta_agm_matrices(acb_theta_newton_matrix(ctx, 0), try, g);
      arf_pos_inf(acb_theta_newton_rho(ctx));
      arf_zero(acb_theta_newton_max(ctx));
      
      for (k = 0; k < n; k++)
	{
	  acb_theta_newton_set_matrix(ctx, k, tau, acb_theta_newton_matrix(ctx, k), prec);
	}
      
      for (k = 0; k < n; k++)
	{
	  agm_radius(rad, acb_theta_newton_mi(ctx, k),
		     acb_theta_newton_nb_bad_steps(ctx, k), prec);
	  
	  /* Propagate radius according to quotients & duplication */
	  propagate_rho(rad, rad, th, acb_theta_newton_matrix(ctx, k), lowprec);
	  arf_min(acb_theta_newton_rho(ctx),
		  acb_theta_newton_rho(ctx), rad, lowprec);
	  
	  /* Update maximum value of Borchardt quotients */
	  arf_div(rad, acb_theta_newton_M0(ctx, k),
		  acb_theta_newton_minf(ctx, k), lowprec, ARF_RND_CEIL);
	  arf_mul_2exp_si(rad, rad, 1);
	  arf_max(acb_theta_newton_max(ctx),
		  acb_theta_newton_max(ctx), rad, lowprec);

	  _arb_vec_clear(mi, nb_bad_steps[k]);
	}

      if (!arf_is_finite(acb_theta_newton_max(ctx))) continue;
      if (!arf_is_positive(acb_theta_newton_rho(ctx))) continue;

      /* Evaluate finite difference */
      acb_theta_cauchy(B2, acb_theta_newton_rho(ctx),
		       acb_theta_newton_max(ctx), 2, lowprec);
      arf_frexp(B2, e, B2);
      exp = fmpz_get_si(e);
      arb_one(eta);
      arb_mul_2exp_si(eta, eta, FLINT_MIN(-exp-n_clog(n,2), -prec/2));
      acb_theta_newton_fd(fd, th, eta, ctx, prec);
      res = acb_mat_inv(fdinv, fd);
      if (!res) continue;
      
      /* Is ||FD^-1||*n*B2*eta less than 1? */
      acb_mat_ninf(norm, fdinv, lowprec);
      arb_mul(bound, norm, B2, lowprec);
      arb_mul_si(bound, bound, n, lowprec);
      arb_mul(bound, bound, eta, lowprec);
      arb_sub_si(test, prod, 1, lowprec);
      if (!arb_is_negative(test)) continue;

      /* Get inv_der */
      arb_mul(bound, bound, norm, lowprec);
      arb_add(bound, bound, norm, lowprec);
      arb_get_ubound_arf(acb_theta_newton_inv_der(ctx), bound, lowprec);
      stop = 1;
    }

  /* Clear */
  acb_mat_clear(half);
  _acb_vec_clear(th, 1<<g);
  arf_clear(rad);
  arf_clear(B2);
  fmpz_clear(e);
  arb_clear(eta);
  acb_mat_clear(fd);
  acb_mat_clear(fdinv);
  arb_clear(norm);
  arb_clear(bound);
  arb_clear(test);
}
