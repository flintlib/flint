
#include "acb_theta.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("agm_sqrt_lowprec....");
  fflush(stdout);
  
  flint_randinit(state);

  /* Test: value of square root should agree; precision remains high */
  for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
      acb_t rt;
      acb_t x;
      arb_t err;
      acb_t rt_low;
      acb_t test;
      
      slong prec = 100 + n_randint(state, 1000);
      slong mag_bits = n_randint(state, 4);
      slong lowprec = ACB_THETA_AGM_LOWPREC;

      acb_init(rt);
      acb_init(x);
      arb_init(err);
      acb_init(rt_low);
      acb_init(test);
      
      acb_randtest_precise(rt, state, prec, mag_bits);
      acb_sqr(x, rt, prec);
      arb_one(err);
      arb_mul_2exp_si(err, err, -lowprec);
      acb_set(rt_low, rt);
      acb_add_error_arb(rt_low, err);

      acb_theta_agm_sqrt_lowprec(test, x, rt_low, prec);

      if (!acb_overlaps(test, rt))
	{
	  flint_printf("FAIL (value)\n");
	  fflush(stdout);
	  flint_abort();
	}

      acb_get_mid(x, test);
      acb_sub(test, test, x, prec);
      acb_abs(err, test, prec);
      arb_mul_2exp_si(err, err, prec - n_pow(2,mag_bits) - ACB_THETA_AGM_GUARD);
      arb_add_si(err, err, -1, prec);

      if (!arb_is_negative(err))
	{
	  flint_printf("FAIL (precision)\n");
	  fflush(stdout);
	  flint_abort();
	}

      acb_clear(rt);
      acb_clear(x);
      arb_clear(err);
      acb_clear(rt_low);
      acb_clear(test);      
    }
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
