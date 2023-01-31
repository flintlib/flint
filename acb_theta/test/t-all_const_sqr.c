
#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("all_const_sqr....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 5 * arb_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 2);
        slong n = 1 << (2 * g);
        slong prec = 2000 + n_randint(state, 2000);
        slong mag_bits = 1 + n_randint(state, 3);

        acb_mat_t tau;
        acb_ptr th2;
        acb_ptr test;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        th2 = _acb_vec_init(n);
        test = _acb_vec_init(n);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);

        /* Force unbalancedness */
        for (j = 0; j < g; j++)
        {
            for (k = 0; k < g; k++)
            {
                acb_mul_2exp_si(acb_mat_entry(tau, j, k),
                                acb_mat_entry(tau, j, k), j + k);
            }
        }

        acb_theta_all_const_sqr(th2, tau, prec);
        acb_theta_naive_all_const(test, tau, prec);
        for (k = 0; k < n; k++)
            acb_sqr(&test[k], &test[k], prec);

        res = 1;
        for (k = 0; k < n; k++)
        {
            if (!acb_overlaps(&th2[k], &test[k]))
                res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(th2, n);
        _acb_vec_clear(test, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
