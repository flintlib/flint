#include "gr.h"

int main(void)
{
    gr_ctx_t G;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong n;

    flint_printf("perm....");
    fflush(stdout);

    for (n = 0; n <= 20; n++)
    {
        gr_ctx_init_perm(G, n);
        gr_test_multiplicative_group(G, 100, flags);
        gr_ctx_clear(G);
    }

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
