
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;
  
    flint_printf("bound....");
    fflush(stdout);
  
    flint_randinit(state);

    /* Test: value of theta should be less than bound */
    for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        slong rad_exp = 1;
        slong n = 1 << (2*g);
        acb_mat_t tau;
        acb_ptr z;
        arf_t rad;
        arf_t bound;
        acb_ptr th;
        arb_t abs;
        arb_t cmp;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        arf_init(rad);
        arf_init(bound);
        th = _acb_vec_init(n);
        arb_init(abs);
        arb_init(cmp);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, rad_exp);
        for (k = 0; k < g; k++)
        {
            acb_randtest_disk(&z[k], &z[k], rad, state, prec);
        }
      
        acb_theta_bound(rad, bound, z, tau, prec);

        if (arf_cmp_si(rad,0) <= 0 || !arf_is_finite(bound))
	{
            flint_printf("Warning: not finite\n");
	}
        else
	{
            for (j = 0; j < g; j++)
	    {
                for (k = 0; k < g; k++)
		{
                    acb_randtest_disk(acb_mat_entry(tau, j, k),
                            acb_mat_entry(tau, j, k), rad, state, prec);
		}
	    }
            for (k = 0; k < g; k++)
	    {
                acb_randtest_disk(&z[k], &z[k], rad, state, prec);
	    }
	}
        acb_theta_naive_all(th, z, tau, prec);

        arb_set_arf(cmp, bound);
        res = 1;
        for (k = 0; k < n; k++)
	{
            acb_abs(abs, &th[k], prec);
            if (arb_gt(abs, cmp)) res = 0;
	}

        if (!res)
	{
            flint_printf("FAIL: theta value is too large\n");
            flint_printf("g = %wd, prec = %wd, tau, z in disk:\n", g, prec);
            acb_mat_printd(tau, 10);
            for (k = 0; k < g; k++)
	    {
                acb_printd(&z[k], 10); flint_printf("\n");
	    }
            flint_printf("rad: "); arf_printd(rad, 10); flint_printf("\n");
            flint_printf("bound: "); arf_printd(bound, 10); flint_printf("\n");
            flint_printf("theta:\n");
            for (k = 0; k < n; k++)
	    {
                acb_printd(&th[k], 10); flint_printf("\n");
	    }	  
            fflush(stdout);
            flint_abort();
	}

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        arf_clear(rad);
        arf_clear(bound);
        _acb_vec_clear(th, n);
        arb_clear(abs);
        arb_clear(cmp);
    }
  
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


      
      
