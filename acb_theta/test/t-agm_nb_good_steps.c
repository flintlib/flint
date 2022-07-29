
#include "acb_theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("agm_nb_good_steps....");
  fflush(stdout);
  
  flint_randinit(state);

  /* Test: for good input, after this number of steps, relative error
     must at most 1/5*2^(-prec) */
  for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
      slong g = 1 + n_randint(state, 4);
      slong prec = ACB_THETA_AGM_BASEPREC + n_randint(state, 1000);
      slong test_prec = prec / (1 + n_randint(state, 10));
      slong n = 1<<g;
      acb_ptr a;
      slong nb_good;
      arf_t rad;
      arb_t err, cmp;
      acb_t diff;
      slong k;
      int res;
      
      a = _acb_vec_init(n);
      arf_init(rad);
      arb_init(err);
      arb_init(cmp);
      acb_init(diff);
      
      arb_one(err);
      arb_div_si(err, err, 50, prec);
      arb_get_lbound_arf(rad, err, prec);
      acb_one(diff);
      
      for (k = 0; k < n; k++) acb_randtest_disk(&a[k], diff, rad, state, prec);
      nb_good = acb_theta_agm_nb_good_steps(rad, g, test_prec);

      /*
      flint_printf("Initial values:\n");
      for (k = 0; k < n; k++)
	{
	  acb_printd(&a[k], 10); flint_printf("\n");
	}
      flint_printf("g = %wd, test_prec = %wd, nb_good = %wd\n", g, test_prec, nb_good);

      if (arf_cmp_2exp_si(rad, -test_prec) > 0)
	{
	  flint_printf("FAIL (error bound)\n");
	  fflush(stdout);
	  flint_abort();
	}
      */

      for (k = 0; k < nb_good; k++) acb_theta_agm_step_good(a, a, g, prec);
      acb_abs(cmp, &a[0], prec);
      arb_mul_arf(cmp, cmp, rad, prec);
      arb_div_si(cmp, cmp, 5, prec);

      res = 1;
      for (k = 1; k < n; k++)
	{
	  acb_sub(diff, &a[k], &a[0], prec);
	  acb_abs(err, diff, prec);
	  if (arb_gt(err, cmp))
	    {
	      flint_printf("Index %wd, err, cmp:\n");
	      arb_printd(err, 10); flint_printf("\n");
	      arb_printd(cmp, 10); flint_printf("\n");
	      res = 0;
	    }
	}

      if (!res)
	{
	  flint_printf("FAIL (values)\n");
	  flint_printf("g = %wd, test_prec = %wd, nb_good = %wd\n", g, test_prec, nb_good);
	  for (k = 0; k < n; k++)
	    {
	      acb_printd(&a[k], 10); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(a, n);
      arf_clear(rad);
      arb_clear(err);
      arb_clear(cmp);
      acb_clear(diff);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
  
