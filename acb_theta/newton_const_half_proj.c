
#include "acb_theta.h"

static void get_der_naive(acb_ptr dth, const acb_mat_t jet_0, const acb_mat_t jet2, slong prec)
{
  /* Depends on what we get from Newton exactly */
}

static void newton_run(acb_ptr th, acb_ptr dth, const acb_mat_t tau,
		       const acb_theta_newton_t ctx, slong prec)
{
  half;
  p;
  current;
  eta;
  nextp;
  fd;
  fdinv;
  res;
  f;
  z0;
  h;
  
  /* Replace tau by its middle point */
  p = start_precision(ctx);
  acb_theta_naive_const_proj(current, half, p);
  remove_errors(current); /* Get new precision */
  nextp = next_precision(ctx, p);
  acb_theta_newton_fd(fd, th, eta, ctx, nextp); /* This is an approximation */
  /* Here we assume that automatic error propagation is better than
     the bound we compute on the derivative; should check this */
  /* Or better: add margin, and double if necessary? */
  res = acb_mat_inv(fdinv, fd, nextp); /* Check that prec remains within bounds */
  /* Assert res */
  /* Should get f(x) as well as fd */
  correct_feedback(z0, tau, nextp);
  newton_iteration(h, fdinv, z0-f, nextp); /* Check that prec remains within bounds */
  _acb_vec_add(current, current, h); /* Here prec depends on size of current */
  /* Get midpoint, update precision, should be larger than p */

  /* At last step, we know fdinv up to some explicit precision, and we
     know it is a close approximation of the true derivative (wrt
     feedback). Can deduce derivatives wrt entries of tau, by
     derivating feedback wrt tau */

  /* Add propagation error for error around tau */
  
}

void acb_theta_newton_const_half_proj(acb_ptr th, acb_ptr dth, const acb_mat_t tau, slong prec)
{
  acb_theta_newton_t ctx;
  acb_mat_t half;
  acb_mat_struct* jet;
  
  slong g = acb_mat_nrows(tau);
  slong n = 1<<g;
  slong lowprec = ACB_THETA_NEWTON_LOWPREC;
  slong baseprec = ACB_THETA_NEWTON_BASEPREC;
  slong k;
  int stop = 0;
  int naive = 0;
  
  acb_theta_newton_init(ctx, g, n);
  acb_mat_init(half, g, g);
  jet = flint_malloc(3 * sizeof(acb_mat_struct));
  acb_mat_init(&jet[0], n, acb_theta_nb_partials(0, g));
  acb_mat_init(&jet[1], n, acb_theta_nb_partials(1, g));
  acb_mat_init(&jet[2], n, acb_theta_nb_partials(2, g));  

  acb_mat_scalar_mul_2exp_si(half, tau, -1);

  /* Attempt to set up newton context */
  while (!stop)
    {
      acb_theta_newton_set_all(ctx, tau, baseprec);
      if (!acb_theta_newton_is_valid_ctx(ctx))
	{	  
	  baseprec *=2;
	  if (baseprec > prec / ACB_THETA_NEWTON_BASEPREC_MAXQ)
	    {
	      stop = 1;
	      naive = 1;
	    }
	}
      else stop = 1;
    }

  if (naive)
    {      
      acb_theta_naive_const_jet(jet, half, 2, prec);
      /* Recover th, dth from this data */
      acb_one(&th[0]);
      for (k = 1; k < n; k++)
	{
	  acb_div(&th[k], acb_mat_entry(&jet[0], k, 0),
		  acb_mat_entry(&jet[0], 0, 0), prec);
	}
      get_der_naive(dth, &jet[0], &jet[2], prec);      
    }

  else /* run Newton */
    {
      newton_run(th, dth, tau, ctx, prec);
    }

  for (k = 0; k < 2; k++) acb_mat_clear(&jet[k]);
  flint_free(jet);
  acb_mat_clear(half);
  acb_theta_newton_clear(ctx);
}
