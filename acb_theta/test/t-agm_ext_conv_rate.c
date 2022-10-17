
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
        arf_t rad, m, M;
        arf_t c1, c2, r;
        arb_t abs;
        arb_t eps;
        acb_t mu;
        slong k, j;

        a = _acb_vec_init(2*n);
        b = _acb_vec_init(2*n);
        acb_init(x);
        arf_init(rad);
        arf_init(m);
        arf_init(M);
        arf_init(c1);
        arf_init(c2);
        arf_init(r);
        arb_init(abs);
        arb_init(eps);
        acb_init(mu);

        /* Generate starting values */
        
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

        /* Get conv rates */
        acb_theta_agm_max_abs(abs, a, n, prec);
        arb_get_ubound_arf(M, abs, prec);
        acb_theta_agm_min_abs(abs, a, n, prec);
        arb_get_lbound_arf(m, abs, prec);
        acb_theta_agm_rel_dist(abs, a+n, n, prec, prec);
        arb_get_ubound_arf(rad, abs, prec);

        acb_theta_agm_ext_conv_rate(c1, c2, r, rad, m, M, prec);
        acb_theta_agm(mu, a+n, NULL, 0, g, prec);

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
            arb_set_arf(abs, r);
            arb_pow_ui(abs, abs, 1<<(j-1), prec);
            arb_mul_arf(abs, abs, c2, prec);      
	    
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
        arf_clear(m);
        arf_clear(M);
        arf_clear(c1);
        arf_clear(c2);
        arf_clear(r);
        arb_clear(abs);
        arb_clear(eps);
        acb_clear(mu);
    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
              

      
              
