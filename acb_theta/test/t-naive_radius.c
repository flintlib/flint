
#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: value of naive_tail should be less than 2^(-prec) */
    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong p = n_randint(state, 10);
        slong prec = 200 + n_randint(state, 1000);
        slong int_prec = prec / (1 + n_randint(state, 10));
        slong lowprec = ACB_THETA_ELD_DEFAULT_PREC;
        slong mag_bits = n_randint(state, 4);
        arb_mat_t Y;
        arf_t eps;
        arf_t R;
        arf_t test;

        arb_mat_init(Y, g, g);
        arf_init(eps);
        arf_init(R);
        arf_init(test);

        arb_mat_randtest_cho(Y, state, prec, mag_bits);
        arf_one(eps);
        arf_mul_2exp_si(eps, eps, -int_prec);
        acb_theta_naive_radius(R, Y, p, eps, lowprec);
        acb_theta_naive_tail(test, R, Y, p, lowprec);

        if (arf_cmp(test, eps) > 0)
        {
            flint_printf("FAIL\n");
            arb_mat_printd(Y, 10);
            flint_printf("g = %wd, p = %wd\n", g, p);
            flint_printf("Objective, tail bound, radius:\n");
            arf_printd(eps, 10);
            flint_printf("\n");
            arf_printd(test, 10);
            flint_printf("\n");
            arf_printd(R, 10);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        arb_mat_clear(Y);
        arf_clear(eps);
        arf_clear(R);
        arf_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
