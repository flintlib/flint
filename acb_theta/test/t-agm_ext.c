
#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("agm_ext....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: overlap at different precisions */
    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 500 + n_randint(state, 1000);
        slong mag_bits = 1 + n_randint(state, 4);
        slong test_prec = prec / (2 + n_randint(state, 9));
        slong n = 1 << g;
        acb_ptr a;
        arf_t rad;
        acb_t x, y, z, t;
        slong k;

        a = _acb_vec_init(2 * n);
        arf_init(rad);
        acb_init(x);
        acb_init(y);
        acb_init(z);
        acb_init(t);

        arf_one(rad);
        arf_mul_2exp_si(rad, rad, -4 - n_randint(state, 4));

        acb_one(x);
        for (k = 0; k < n; k++)
            acb_randtest_disk(&a[k], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a, a, n, x, prec);

        arf_one(rad);
        arf_mul_2exp_si(rad, rad, -4 - n_randint(state, 4));

        acb_one(x);
        for (k = 0; k < n; k++)
            acb_randtest_disk(&a[k + n], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a + n, a + n, n, x, prec);

        acb_theta_agm_ext(x, y, a, NULL, 0, g, test_prec);
        acb_theta_agm_ext(z, t, a, NULL, 0, g, prec);

        if (!acb_overlaps(x, z) || !acb_overlaps(y, t))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, test_prec = %wd\n", g, test_prec);
            for (k = 0; k < 2 * n; k++)
            {
                acb_printd(&a[k], 10);
                flint_printf("\n");
            }
            flint_printf("agms:\n");
            acb_printd(x, 10);
            flint_printf("\n");
            acb_printd(y, 10);
            flint_printf("\n");
            acb_printd(z, 10);
            flint_printf("\n");
            acb_printd(t, 10);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _acb_vec_clear(a, 2 * n);
        arf_clear(rad);
        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
