
#include "acb_theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("agm_ext_step....");
  fflush(stdout);
  
  flint_randinit(state);

  /* Test: should coincide with theta duplication */
  for (iter = 0; iter < 20 * arb_test_multiplier(); iter++)
    {
      slong g = 1 + n_randint(state, 3);
      slong prec = 100 + n_randint(state, 500);
      slong mag_bits = n_randint(state, 2);
      slong rad_exp = -5;
      slong n = 1<<g;
      slong lowprec = ACB_THETA_AGM_LOWPREC;
      
      acb_mat_t tau;
      acb_ptr z;
      arf_t rad;
      acb_ptr th;
      acb_ptr th_sqr;
      acb_ptr th_dupl;
      acb_ptr test;
      arb_t err;
      slong k;
      int pos;
      int res;

      acb_mat_init(tau, g, g);
      arf_init(rad);
      z = _acb_vec_init(2*g);
      th = _acb_vec_init(2*n);
      th_sqr = _acb_vec_init(2*n);
      th_dupl = _acb_vec_init(2*n);
      test = _acb_vec_init(2*n);

      acb_siegel_randtest(tau, state, prec, mag_bits);
      arf_one(rad);
      arf_mul_2exp_si(rad, rad, rad_exp);
      for (k = 0; k < g; k++) acb_randtest_disk(&z[k], &z[k], rad, state, prec);

      acb_theta_naive(th, z, 2, tau, prec);

      /*
      flint_printf("g = %wd, prec = %wd, tau, z:\n", g, prec);
      acb_mat_printd(tau, 10);
      for (k = 0; k < g; k++)
	{
	  acb_printd(&z[k], 10); flint_printf("\n");
	}
      flint_printf("theta:\n");
      for (k = 0; k < 2*n; k++)
	{
	  acb_printd(&th[k], 10); flint_printf("\n");
	}
      */      
      
      for (k = 0; k < 2*n; k++) acb_sqr(&th_sqr[k], &th[k], prec);
      arb_one(err);
      arb_mul_2exp_si(err, err, -lowprec);
      for (k = 0; k < 2*n; k++) acb_add_error_arb(&th[k], err);
      
      acb_mat_scalar_mul_2exp_si(tau, tau, 1);
      
      acb_theta_naive(th_dupl, z, 2, tau, prec);
      for (k = 0; k < 2*n; k++) acb_sqr(&th_dupl[k], &th_dupl[k], prec);

      pos = 1;
      for (k = 0; k < 2*n; k++)
	{
	  if (!arb_is_positive(acb_realref(&th[k])))
	    {
	      pos = 0;
	      break;
	    }
	}

      if (pos && (iter%2 == 0)) acb_theta_agm_ext_step_good(test, th_sqr, g, prec);
      else acb_theta_agm_ext_step_bad(test, th_sqr, th, g, prec);
      
      res = 1;
      for (k = 0; k < n; k++)
	{
	  if (!acb_overlaps(&test[k], &th_dupl[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (overlap)\n");	  
	  flint_printf("g = %wd, prec = %wd, tau, z:\n", g, prec);
	  acb_mat_printd(tau, 10);
	  for (k = 0; k < g; k++)
	    {
	      acb_printd(&z[k], 10); flint_printf("\n");
	    }
	  flint_printf("theta:\n");
	  for (k = 0; k < 2*n; k++)
	    {
	      acb_printd(&th[k], 10); flint_printf("\n");
	    }
	  flint_printf("dupl:\n");
	  for (k = 0; k < 2*n; k++)
	    {
	      acb_printd(&th_dupl[k], 10); flint_printf("\n");
	    }
	  flint_printf("test:\n");
	  for (k = 0; k < 2*n; k++)
	    {
	      acb_printd(&test[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      acb_mat_clear(tau);
      arf_clear(rad);
      _acb_vec_clear(z, 2*g);
      _acb_vec_clear(th, 2*n);
      _acb_vec_clear(th_sqr, 2*n);
      _acb_vec_clear(th_dupl, 2*n);
      _acb_vec_clear(test, 2*n);
    }  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}  
