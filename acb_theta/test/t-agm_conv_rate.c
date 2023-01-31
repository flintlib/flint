
#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_conv_rate....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: convergence rate works for first few steps */
    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 500 + n_randint(state, 1000);
        slong mag_bits = 1 + n_randint(state, 4);
        slong n = 1 << g;
        acb_ptr a;
        acb_t x;
        arf_t rad;
        arf_t c, r;
        arb_t abs;
        arb_t eps;
        arf_t err;
        slong k, j;

        a = _acb_vec_init(n);
        acb_init(x);
        arf_init(rad);
        arf_init(c);
        arf_init(r);
        arb_init(abs);
        arb_init(eps);
        arf_init(err);

        arb_one(eps);
        arb_div_si(eps, eps, 16 + n_randint(state, 100), prec);
        arb_get_lbound_arf(rad, eps, prec);

        acb_one(x);
        for (k = 0; k < n; k++)
            acb_randtest_disk(&a[k], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a, a, n, x, prec);

        acb_theta_agm_step_good(a, a, g, prec);
        acb_theta_agm_rel_dist(eps, a, n, prec, prec);
        arb_get_lbound_arf(rad, eps, prec);
        acb_theta_agm_conv_rate(c, r, rad, prec);
        for (j = 0; j < 5; j++)
        {
            /* Get current relative error */
            arb_zero(eps);
            for (k = 0; k < n; k++)
            {
                acb_sub(x, &a[0], &a[k], prec);
                acb_abs(abs, x, prec);
                arb_max(eps, eps, abs, prec);
            }
            acb_abs(abs, &a[0], prec);
            arb_div(eps, eps, abs, prec);

            /* Get predicted error */
            arb_set_arf(abs, r);
            arb_pow_ui(abs, abs, 1 << j, prec);
            arb_mul_arf(abs, abs, c, prec);

            if (arb_lt(abs, eps))
            {
                flint_printf("FAIL (error bound)\n");
                flint_printf("At step %wd, predicted and real error \n", j);
                arb_printd(abs, 10);
                flint_printf("\n");
                arb_printd(eps, 10);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            acb_theta_agm_step_good(a, a, g, prec);
        }

        _acb_vec_clear(a, n);
        acb_clear(x);
        arf_clear(rad);
        arf_clear(c);
        arf_clear(r);
        arb_clear(abs);
        arb_clear(eps);
        arf_clear(err);
    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
