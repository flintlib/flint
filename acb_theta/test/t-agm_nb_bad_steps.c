
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;
  
    flint_printf("agm_nb_bad_steps....");
    fflush(stdout);
  
    flint_randinit(state);

    /* Test: after scalar multiplication, theta values must be close to 1 */
    for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_AGM_BASEPREC + n_randint(state, 1000);
        slong mag_bits = n_randint(state, 2);
        slong n = 1 << g;
        acb_mat_t tau;
        slong nb_bad;
        acb_ptr th;
        acb_t diff;
        arb_t err;
        arb_t cmp;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        th = _acb_vec_init(n);
        acb_init(diff);
        arb_init(err);
        arb_init(cmp);
        
        acb_siegel_randtest(tau, state, prec, mag_bits);
        nb_bad = acb_theta_agm_nb_bad_steps(tau, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, nb_bad);
        acb_theta_naive_const(th, tau, prec);

        /* Theta values must be relatively close to &th[0], at 1/20 */
        acb_abs(cmp, &th[0], prec);
        arb_div_si(cmp, cmp, 20, prec);
        res = 1;
        for (k = 1; k < n; k++)
	{
            acb_sub(diff, &th[k], &th[0], prec);
            acb_abs(err, diff, prec);
            if (arb_gt(err, cmp))
	    {
                res = 0;
                break;
	    }
	}
        
        if (!res)
	{
	  flint_printf("FAIL (theta values are not close)\n");
          flint_printf("tau:\n");
          acb_mat_printd(tau, 10);
          flint_printf("theta values:\n");
          for (k = 0; k < n; k++)
          {
              acb_printd(&th[k], 10); flint_printf("\n");
          }
	  fflush(stdout);
	  flint_abort();
	}
        
        acb_mat_clear(tau);
        _acb_vec_clear(th, n);
        acb_clear(diff);
        arb_clear(err);
        arb_clear(cmp);      
    }
  
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
