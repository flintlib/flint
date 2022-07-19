
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("const_naive....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with const_ind_naive; duplication formula */
    for (iter = 0; iter < 200 * arb_test_multiplier(); iter++)
      {	
	slong g = 1 + n_randint(state, 4);
	slong nb = n_pow(2,g);
	acb_mat_t tau;
	acb_ptr th;
	acb_ptr th_test;
	ulong ab;
	slong prec = 20 + n_randint(state, 2000);
	slong mag_bits = n_randint(state, 2);
	int res;
	slong k;
	
	acb_mat_init(tau, g, g);
	th = _acb_vec_init(nb);
	th_test = _acb_vec_init(nb);

	acb_siegel_randtest(tau, state, prec, mag_bits);
	
	flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
	acb_mat_printd(tau, 10);
	
	for (ab = 0; ab < nb; ab++)
	  {
	    acb_theta_const_ind_naive(&th_test[ab], ab, tau, prec);
	  }
	acb_theta_const_naive(th, tau, prec);

	res = 1;
	for (k = 0; k < nb; k++)
	  {
	    if (!acb_overlaps(&th[k], &th_test[k])) res = 0;
	  }
	if (!res)
	  {
	    flint_printf("FAIL: overlap\n");
	    flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
	    acb_mat_printd(tau, 10);
	    flint_printf("th[k], th_test[k]:\n");
	    for (k = 0; k < nb; k++)
	      {
		acb_printd(&th[k], 10); flint_printf("\n");
		acb_printd(&th_test[k], 10); flint_printf("\n");
	      }
	    fflush(stdout);
	    flint_abort();
	  }
	
	acb_theta_duplication(th_test, th, g, prec);
	acb_mat_scalar_mul_2exp_si(tau, tau, 1);	
	acb_theta_const_naive(th, tau, prec);
	for (k = 0; k < nb; k++) acb_sqr(&th[k], &th[k], prec);
	
	res = 1;
	for (k = 0; k < nb; k++)
	  {
	    if (!acb_overlaps(&th[k], &th_test[k])) res = 0;
	  }
	if (!res)
	  {
	    flint_printf("FAIL: duplication\n");
	    flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
	    acb_mat_printd(tau, 10);
	    flint_printf("th[k], th_test[k]:\n");
	    for (k = 0; k < nb; k++)
	      {
		acb_printd(&th[k], 10); flint_printf("\n");
		acb_printd(&th_test[k], 10); flint_printf("\n");
	      }
	    fflush(stdout);
	    flint_abort();
	  }

	acb_mat_clear(tau);
	_acb_vec_clear(th, nb);
	_acb_vec_clear(th_test, nb);	
      }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
