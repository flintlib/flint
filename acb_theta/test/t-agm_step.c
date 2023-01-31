
#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_step....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: should coincide with theta duplication */
    for (iter = 0; iter < 20 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        slong n = 1 << g;
        slong lowprec = ACB_THETA_AGM_LOWPREC;

        acb_mat_t tau;
        acb_ptr th;
        acb_ptr th_sqr;
        acb_ptr th_dupl;
        acb_ptr test;
        arb_t err;
        slong k;
        int pos;
        int res;

        acb_mat_init(tau, g, g);
        th = _acb_vec_init(n);
        th_sqr = _acb_vec_init(n);
        th_dupl = _acb_vec_init(n);
        test = _acb_vec_init(n);
        arb_init(err);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_naive_const(th, tau, prec);

        /*
           flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
           acb_mat_printd(tau, 10);
           flint_printf("theta:\n");
           for (k = 0; k < n; k++)
           {
           acb_printd(&th[k], 10); flint_printf("\n");
           }
         */

        for (k = 0; k < n; k++)
            acb_sqr(&th_sqr[k], &th[k], prec);
        arb_one(err);
        arb_mul_2exp_si(err, err, -lowprec);
        for (k = 0; k < n; k++)
            acb_add_error_arb(&th[k], err);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        acb_theta_naive_const(th_dupl, tau, prec);
        for (k = 0; k < n; k++)
            acb_sqr(&th_dupl[k], &th_dupl[k], prec);

        pos = 1;
        for (k = 0; k < n; k++)
        {
            if (!arb_is_positive(acb_realref(&th[k])))
            {
                pos = 0;
                break;
            }
        }

        if (pos && iter % 2 == 0)
            acb_theta_agm_step_good(test, th_sqr, g, prec);
        else
            acb_theta_agm_step_bad(test, th_sqr, th, g, prec);

        res = 1;
        for (k = 0; k < n; k++)
        {
            if (!acb_overlaps(&test[k], &th_dupl[k]))
                res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 10);
            flint_printf("theta:\n");
            for (k = 0; k < n; k++)
            {
                acb_printd(&th[k], 10);
                flint_printf("\n");
            }
            flint_printf("dupl:\n");
            for (k = 0; k < n; k++)
            {
                acb_printd(&th_dupl[k], 10);
                flint_printf("\n");
            }
            flint_printf("test:\n");
            for (k = 0; k < n; k++)
            {
                acb_printd(&test[k], 10);
                flint_printf("\n");
            }
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(th, n);
        _acb_vec_clear(th_sqr, n);
        _acb_vec_clear(th_dupl, n);
        _acb_vec_clear(test, n);
        arb_clear(err);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
