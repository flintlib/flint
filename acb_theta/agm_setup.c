
#include "acb_theta.h"

void acb_theta_agm_setup(fmpz_mat_struct* Ni, slong* nb_bad_steps, acb_ptr all_roots,
			 arb_t rho, arb_t M, arb_t Binv,
			 const acb_mat_t tau, slong prec)
{
  acb_mat_t half;
  acb_ptr th;
  slong index;
  acb_ptr all_r0;
  arb_t M0, minf;
  arb_ptr mi;
  arb_t term, prod, rad, norm;
  arb_t B2;
  fmpz_t e;
  slong exp;
  arb_t eta;
  acb_mat_t fd, fdinv;
  
  slong lowprec = ACB_THETA_AGM_LOWPREC;  
  slong g = acb_mat_nrows(tau);
  slong n = (1 << g) - 1;
  slong k, j;
  int stop = 0;
  int res;
  int try = 0;

  /* Init */

  acb_mat_scalar_mul_2exp_si(half, tau, -1);
  acb_theta_naive_const(th, half, prec);

  while (!stop && (try < ACB_THETA_AGM_NB_MATRIX_SETUPS))
    {
      try++;
      acb_theta_agm_matrices(Ni, try, g);
      arb_pos_inf(rho);
      arb_zero(M);
      
      for (k = 0; k < n; k++)
	{
	  nb_bad_steps[k] = acb_theta_agm_nb_bad_steps(tau, &Ni[k], prec);
	  nb_total += nb_bad_steps[k];
	}
      
      /* This won't work, make structure and update size... */
      index = 0;
      for (k = 0; k < n; k++)
	{
	  mi = _arb_vec_init(nb_bad_steps[k]);

	  /* Collect info on AGM sequence */
	  acb_theta_agm_collect(all_roots + index, M0, minf, mi,
				nb_bad_steps[k], tau, N, lowprec);
	  index += nb_bad_steps[k];

	  /* Compute radius */
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

	  /* Propagate radius according to quotients & duplication */
	  acb_theta_propagate_radius(rad, rad, th, &Ni[k], lowprec);
	  arb_min(rho, rad, lowprec);
	  
	  /* Update maximum value of Borchardt quotients */
	  arb_div(term, M0, minf, lowprec);
	  arb_mul_2exp_si(term, term, 1);
	  arb_max(M, M, term, lowprec);

	  _arb_vec_clear(mi, nb_bad_steps[k]);
	}

      /* Now we know everything except Binv. Evaluate finite difference */
      acb_theta_cauchy(B2, rho, M, 2, lowprec);
      arf_frexp(B2, e, B2);
      exp = fmpz_get_si(e);
      arb_one(eta);
      arb_mul_2exp_si(eta, eta, FLINT_MIN(-exp-n_clog(n,2), -prec/2));

      /* Is FD invertible? */
      acb_theta_agm_fd(fd, th, Ni, eta, prec);
      res = acb_mat_inv(fdinv, fd);
      if (!res) continue;
      
      /* Is ||FD^-1||*n*B2*eta less than 1? */
      acb_mat_norm(norm, fdinv, lowprec);
      arb_mul(prod, norm, B2, lowprec);
      arb_mul_si(prod, prod, n, lowprec);
      arb_mul(prod, prod, eta, lowprec);
      arb_sub_si(term, prod, 1, lowprec);
      if (!arb_is_negative(term)) continue;

      /* Get Binv */
      arb_mul(Binv, prod, norm, lowprec);
      arb_add(Binv, Binv, norm, lowprec);
      stop = 1;
    }

  /* Clear */
}
