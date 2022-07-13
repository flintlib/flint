#include "gr.h"

int main()
{
    gr_ctx_t ZZn, ZZnx;
    slong iter;
    int flags = 0;
    flint_rand_t state;

    flint_printf("mpoly_nmod8....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 50; iter++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);

        gr_ctx_init_mpoly(ZZnx, ZZn, n_randint(state, 3), mpoly_ordering_randtest(state));
        ZZnx->size_limit = 100;

        gr_test_ring(ZZnx, 100, flags);
        gr_ctx_clear(ZZnx);

        gr_ctx_clear(ZZn);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
