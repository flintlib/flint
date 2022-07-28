
#include "acb_theta.h"

void acb_theta_newton_const_half_proj(acb_ptr th, const acb_mat_t tau, slong prec)
{
  acb_theta_agm_ctx_t ctx;
  acb_mat_t half;
  /* acb_mat_struct* jet; */
  
  slong g = acb_mat_nrows(tau);
  slong n = 1<<g;
  slong lowprec = ACB_THETA_AGM_CTX_LOWPREC;
  slong baseprec = ACB_THETA_AGM_CTX_BASEPREC;
  slong k;
  int stop = 0;
  int naive = 0;
  
  acb_theta_agm_ctx_init(ctx, g, n);
  acb_mat_init(half, g, g);

  /*
  jet = flint_malloc(3 * sizeof(acb_mat_struct));
  acb_mat_init(&jet[0], n, acb_theta_nb_partials(0, g));
  acb_mat_init(&jet[1], n, acb_theta_nb_partials(1, g));
  acb_mat_init(&jet[2], n, acb_theta_nb_partials(2, g));  
  */

  acb_mat_scalar_mul_2exp_si(half, tau, -1);

  /* Attempt to set up newton context */
  while (!stop)
    {
      acb_theta_agm_ctx_set_all(ctx, tau, baseprec);
      if (!acb_theta_agm_ctx_is_valid(ctx))
	{	  
	  baseprec *=2;
	  if (baseprec > prec / ACB_THETA_AGM_CTX_BASEPREC_MAXQ)
	    {
	      stop = 1;
	      naive = 1;
	    }
	}
      else stop = 1;
    }

  if (naive)
    {      
      /*acb_theta_naive_const_jet(jet, half, 2, prec);
      acb_one(&th[0]);
      for (k = 1; k < n; k++)
	{
	  acb_div(&th[k], acb_mat_entry(&jet[0], k, 0),
		  acb_mat_entry(&jet[0], 0, 0), prec);
	}
	get_der_naive(dth, &jet[0], &jet[2], prec);      */
      acb_theta_naive_const_proj(th, half, prec);
    }

  else /* run Newton */
    {
      acb_theta_newton_run(th, tau, ctx, prec);
    }

  /*for (k = 0; k < 2; k++) acb_mat_clear(&jet[k]);
    flint_free(jet);*/
  acb_mat_clear(half);
  acb_theta_newton_clear(ctx);
}
