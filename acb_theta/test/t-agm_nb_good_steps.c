
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;
  
    flint_printf("agm_nb_good_steps....");
    fflush(stdout);
  
    flint_randinit(state);

    /* Test: overlap at different precisions */
    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 500 + n_randint(state, 1000);
        slong mag_bits = 1 + n_randint(state, 4);
        slong test_prec = prec / (2 + n_randint(state, 9));
        slong n = 1<<g;
        acb_ptr a;
        arf_t rad;
        arf_t c, e;
        slong nb_good;
        acb_t x, y;
        slong k;
      
        a = _acb_vec_init(n);
        arf_init(rad);
        arf_init(c);
        arf_init(e);
        acb_init(x);
        acb_init(y);
      
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, -4);
      
        acb_one(x);      
        for (k = 0; k < n; k++) acb_randtest_disk(&a[k], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a, a, n, x, prec);
      
        acb_theta_agm_conv_rate(c, e, a, g, prec);
      
        nb_good = acb_theta_agm_nb_good_steps(c, e, test_prec);
        acb_theta_agm(x, a, NULL, 0, nb_good, g, test_prec);
        nb_good = acb_theta_agm_nb_good_steps(c, e, prec);
        acb_theta_agm(y, a, NULL, 0, nb_good, g, prec);   

        if (!acb_overlaps(x, y))
	{
            flint_printf("FAIL (values)\n");   
            flint_printf("g = %wd, test_prec = %wd, nb_good = %wd\n", g, test_prec, nb_good);
            for (k = 0; k < n; k++)
            {
                acb_printd(&a[k], 10); flint_printf("\n");
            }
            flint_printf("agm:\n");
            acb_printd(x, 10); flint_printf("\n");
            acb_printd(y, 10); flint_printf("\n");
            fflush(stdout);
            flint_abort();
	}

        _acb_vec_clear(a, n);
        arf_clear(rad);
        arf_clear(c);
        arf_clear(e);
        acb_clear(x);
        acb_clear(y);
    }
  
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
      
  
