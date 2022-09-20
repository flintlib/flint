
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("newton_const_sqr....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 2 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = 1<<g;
        acb_mat_t tau;
        acb_ptr th2;
        acb_ptr th2_test;
        slong prec = (2 + n_randint(state, 4-g)) * ACB_THETA_AGM_BASEPREC
            + 100 * n_randint(state, 20);
        int res;
        slong k;

        acb_mat_init(tau, g, g);
        th2 = _acb_vec_init(nb);
        th2_test = _acb_vec_init(nb);
        
        acb_siegel_randtest_fund(tau, state, prec);
        
        flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
        acb_mat_printd(tau, 10); flint_printf("\n");

        flint_printf("Naive...\n");
        acb_theta_naive_const(th2_test, tau, prec);
        flint_printf("Done.\n");
        
        acb_theta_newton_const_sqr(th2, tau, prec);
        
        res = 1;
        for (k = 0; k < nb; k++)
        {
            acb_sqr(&th2_test[k], &th2_test[k], prec);
            if (!acb_overlaps(&th2_test[k], &th2[k])) res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("th[k], th_test[k]:\n");
            for (k = 0; k < nb; k++)
            {
                acb_printd(&th2[k], 10); flint_printf("\n");
                acb_printd(&th2_test[k], 10); flint_printf("\n\n");
            }
            flint_printf("Diff:\n");
            for (k = 0; k < nb; k++)
            {
                acb_sub(&th2[k], &th2[k], &th2_test[k], prec);
                acb_printd(&th2[k], 10); flint_printf("\n");
            }
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(th2, nb);
        _acb_vec_clear(th2_test, nb);        
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

