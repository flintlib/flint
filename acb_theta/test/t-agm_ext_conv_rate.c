
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;
  
    flint_printf("agm_ext_conv_rate....");
    fflush(stdout);
  
    flint_randinit(state);

    /* Test: convergence rate works for first few steps */
    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 500 + n_randint(state, 1000);
        slong mag_bits = 1 + n_randint(state, 4);
        slong n = 1<<g;
        acb_ptr a, b;
        acb_t x;
        arf_t rad;
        arf_t c_ext, c, e;
        arb_t abs;
        arb_t eps;
        acb_t mu;
        arf_t err;
        slong nb_good;
        slong k, j;

        a = _acb_vec_init(2*n);
        b = _acb_vec_init(2*n);
        acb_init(x);
        arf_init(rad);
        arf_init(c_ext);
        arf_init(c);
        arf_init(e);
        arb_init(abs);
        arb_init(eps);
        acb_init(mu);
        arf_init(err);

        arb_one(eps);
        arb_div_si(eps, eps, 16 + n_randint(state, 100), prec);
        arb_get_lbound_arf(rad, eps, prec);

        acb_one(x);
        for (k = 0; k < n; k++) acb_randtest_disk(&a[k], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a, a, n, x, prec);
        
        arb_one(eps);
        arb_div_si(eps, eps, 16 + n_randint(state, 100), prec);
        arb_get_lbound_arf(rad, eps, prec);

        acb_one(x);
        for (k = 0; k < n; k++) acb_randtest_disk(&a[k+n], x, rad, state, prec);
        arb_randtest_pos(acb_realref(x), state, prec, mag_bits);
        _acb_vec_scalar_mul(a+n, a+n, n, x, prec);

        acb_theta_agm_ext_conv_rate(c_ext, c, e, a, g, prec);
        nb_good = acb_theta_agm_nb_good_steps(c, e, prec);
        acb_theta_agm(mu, a+n, NULL, 0, nb_good, g, prec);

        acb_theta_agm_ext_step_good(a, a, g, prec);

        for (j = 1; j < 5; j++)
        {
            /* Get q_(n+1) */
            acb_theta_agm_ext_step_good(b, a, g, prec);
            acb_sqr(x, &b[0], prec);
            acb_div(x, x, &a[0], prec);
            acb_div(x, x, mu, prec);

            acb_sub_si(x, x, 1, prec);
            acb_abs(eps, x, prec);

            /* Get predicted error */
            arb_set_arf(abs, e);
            arb_pow_ui(abs, abs, 1<<(j-1), prec);
            arb_mul_arf(abs, abs, c_ext, prec);
            
            if (arb_lt(abs, eps))
            {
                flint_printf("FAIL (error bound)\n");
                flint_printf("At step %wd, predicted and real error \n", j);
                arb_printd(abs, 10); flint_printf("\n");
                arb_printd(eps, 10); flint_printf("\n");                                
                flint_printf("q = "); acb_printd(x, 10); flint_printf("\n");
                flint_printf("Values:\n");
                acb_printd(&b[0], 10); flint_printf("\n");
                acb_printd(&a[0], 10); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
            
            _acb_vec_set(a, b, 2*n);
        }
      
        _acb_vec_clear(a, 2*n);
        _acb_vec_clear(b, 2*n);
        acb_clear(x);
        arf_clear(rad);
        arf_clear(c_ext);
        arf_clear(c);
        arf_clear(e);
        arb_clear(abs);
        arb_clear(eps);
        acb_clear(mu);
        arf_clear(err);
    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
              

      
              
