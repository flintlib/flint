
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive_const_proj....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: theta is a modular form */
    for (iter = 0; iter < 50 * arb_test_multiplier(); iter++)
    {	
	slong g = 1 + n_randint(state, 2);
	slong nb = n_pow(2,g);
	acb_mat_t tau, Ntau;
	acb_ptr th, th_dupl, th_test;
        acb_t scal;
        fmpz_mat_t N;
	slong prec = 20 + n_randint(state, 300);
	slong mag_bits = n_randint(state, 2);
	int res;
	slong k;
        
	acb_mat_init(tau, g, g);
        acb_mat_init(Ntau, g, g);
	th = _acb_vec_init(nb);
	th_dupl = _acb_vec_init(nb*nb);
	th_test = _acb_vec_init(nb);
        acb_init(scal);
        fmpz_mat_init(N, 2*g, 2*g);

        acb_siegel_randtest_fund(tau, state, prec);
        fmpz_mat_randtest_sp(N, state, mag_bits);
        
        acb_theta_naive_const_proj(th, tau, prec);
        acb_theta_dupl_all_const(th_dupl, th, g, prec);
        acb_theta_transform_sqr_proj(th_test, th_dupl, N, prec);

        acb_inv(scal, &th_test[0], prec);
        _acb_vec_scalar_mul(th_test, th_test, nb, scal, prec);

        acb_mat_scalar_mul_2exp_si(Ntau, tau, 1);
        acb_siegel_transform(Ntau, N, Ntau, prec);
        acb_theta_naive_const_proj(th, Ntau, prec);
        for (k = 0; k < nb; k++) acb_sqr(&th[k], &th[k], prec);
        
        acb_inv(scal, &th[0], prec);
        _acb_vec_scalar_mul(th, th, nb, scal, prec);

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
	    flint_printf("th_test[k], th[k]:\n");
	    for (k = 0; k < nb; k++)
            {
		acb_printd(&th_test[k], 100); flint_printf("\n");
		acb_printd(&th[k], 100); flint_printf("\n");
                flint_printf("\n");
            }
	    fflush(stdout);
	    flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(Ntau);
        _acb_vec_clear(th, nb);
        _acb_vec_clear(th_dupl, nb*nb);
        _acb_vec_clear(th_test, nb);
        acb_clear(scal);
        fmpz_mat_clear(N);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

