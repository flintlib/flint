
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("siegel_reduce....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: mat is symplectic; upper left imag entry is not less than 1/2,
       and real part is reduced */
    for (iter = 0; iter < 200 * arb_test_multiplier(); iter++)
    {
	slong g = 1 + n_randint(state, 4);
	slong prec = 100 + n_randint(state, 500);
	slong mag_bits = 1 + n_randint(state, 5);
	
	acb_mat_t tau;
	acb_mat_t res;
	fmpz_mat_t mat;
	arb_t test;
	int r = 1;

	acb_mat_init(tau, g, g);
	acb_mat_init(res, g, g);
	fmpz_mat_init(mat, 2*g, 2*g);
	arb_init(test);

	acb_siegel_randtest(tau, state, prec, mag_bits);
	acb_siegel_reduce(res, mat, tau, prec);

	if (!fmpz_mat_is_sp(mat))
	{
	    flint_printf("FAIL (symplectic)\n");
	    fmpz_mat_print(mat);
	    fflush(stdout);
	    flint_abort();
	}

	acb_abs(test, acb_mat_entry(res, 0, 0), prec);
	arb_mul_2exp_si(test, test, 1);
	arb_sub_si(test, test, 1, prec);
	if (arb_is_negative(test)) r = 0;

	arb_abs(test, acb_realref(acb_mat_entry(res, 0, 0)));
	arb_sub_si(test, test, 1, prec);
	if (arb_is_positive(test)) r = 0;

	if (!r)
	{
	    flint_printf("FAIL (not reduced):\n");
	    acb_mat_printd(tau, 10);
	    fmpz_mat_print_pretty(mat); flint_printf("\n");
	    acb_mat_printd(res, 10);
	    fflush(stdout);
	    flint_abort();
	}
	
	acb_mat_clear(tau);
	acb_mat_clear(res);
	fmpz_mat_clear(mat);
	arb_clear(test);
	
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
