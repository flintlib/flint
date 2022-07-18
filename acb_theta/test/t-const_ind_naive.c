
#include "acb_theta.h"
#include "acb_modular.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("const_ind_naive....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with genus 1 theta */
    for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
      {
	slong g = 1;
	acb_mat_t tau;
	acb_t t;
	arb_t eps;
	acb_t th;
	acb_t th1, th2, th3, th4, z;
	ulong ab;
	slong prec = 20 + n_randint(state, 100);
	slong mag_bits = n_randint(state, 10);
	
	acb_mat_init(tau, g, g);
	acb_init(t);
	acb_init(th);
	acb_init(th1);
	acb_init(th2);
	acb_init(th3);
	acb_init(th4);
	acb_init(z);
	arb_init(eps);

	arb_one(eps);
	arb_mul_2exp_si(eps, eps, -n_randint(state, 10)); /* Not too small */

	arb_randtest_precise(acb_realref(t), state, prec, mag_bits);
	arb_randtest_precise(acb_imagref(t), state, prec, mag_bits);
	arb_sqr(acb_imagref(t), acb_imagref(t), prec);
	arb_add(acb_imagref(t), acb_imagref(t), eps, prec);

	acb_set(acb_mat_entry(tau, 0, 0), t);

	acb_zero(z);
	acb_modular_theta(th1, th2, th3, th4, z, t, prec);

	flint_printf("tau:\n");
	acb_mat_printd(tau, 10); flint_printf("\n");
	
	flint_printf("Genus 1 theta:\n");
	acb_printd(th1, 30); flint_printf("\n");
	acb_printd(th2, 30); flint_printf("\n");
	acb_printd(th3, 30); flint_printf("\n");
	acb_printd(th4, 30); flint_printf("\n");

	flint_printf("General naive algorithm:\n");
	for (ab = 0; ab < 4; ab++)
	  {
	    acb_theta_const_ind_naive(th, ab, tau, prec);
	    acb_printd(th, 30); flint_printf("\n");	    
	  }
	
	acb_mat_clear(tau);
	acb_clear(t);
	acb_clear(th);
	acb_clear(th1);
	acb_clear(th2);
	acb_clear(th3);
	acb_clear(th4);
	acb_clear(z);
	arb_clear(eps);
      }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
