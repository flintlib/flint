#include "gr.h"

int main()
{
    gr_ctx_t G;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong q;

    flint_printf("dirichlet....");
    fflush(stdout);

    for (q = 1; q <= 100; q++)
    {
        GR_MUST_SUCCEED(gr_ctx_init_dirichlet_group(G, q));
        gr_test_multiplicative_group(G, 100, flags);
        gr_ctx_clear(G);
    }

    {
        flint_rand_t state;
        slong iter;
        flint_randinit(state);

        for (iter = 0; iter < 100; iter++)
        {
            q = n_randtest(state);
            if (gr_ctx_init_dirichlet_group(G, q) == GR_SUCCESS)
            {
                gr_test_multiplicative_group(G, 100, flags);
                gr_ctx_clear(G);
            }
        }

        flint_randclear(state);
    }

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
