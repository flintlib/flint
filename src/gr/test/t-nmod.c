#include "ulong_extras.h"
#include "gr.h"

int main()
{
    gr_ctx_t ZZn;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong n;

    flint_rand_t state;
    flint_randinit(state);

    flint_printf("nmod....");
    fflush(stdout);

    for (n = 1; n < 100; n++)
    {
        gr_ctx_init_nmod(ZZn, n_randtest_not_zero(state));
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    flint_randclear(state);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
