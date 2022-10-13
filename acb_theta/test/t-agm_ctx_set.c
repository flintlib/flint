
#include "acb_theta.h"

int main()
{
    slong iter;
    flint_rand_t state;
  
    flint_printf("agm_ctx_set....");
    fflush(stdout);
  
    flint_randinit(state);

    /* Test: agm context must be valid for g=1, 2 and tau in fundamental domain */
    for (iter = 0; iter < 5 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong prec = ACB_THETA_AGM_BASEPREC + n_randint(state, 1000);
        acb_mat_t tau;
        acb_theta_agm_ctx_t ctx;
        int res;

        acb_mat_init(tau, g, g);

        acb_siegel_randtest_fund(tau, state, prec);
        acb_theta_agm_ctx_init(ctx, tau);
        
        res = acb_theta_agm_ctx_set(ctx, prec);

        if (!res)
	{
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
	}
      
        acb_mat_clear(tau);
        acb_theta_agm_ctx_clear(ctx);
    }
  
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
      
  
