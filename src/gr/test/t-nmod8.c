#include "gr.h"

int main(void)
{
    gr_ctx_t ZZn;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong n;

    flint_printf("nmod8....");
    fflush(stdout);

    for (n = 1; n < 256; n++)
    {
        gr_ctx_init_nmod8(ZZn, n);
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    gr_ctx_init_nmod8(ZZn, 107);
    gr_test_ring(ZZn, 10000, flags);
    gr_ctx_clear(ZZn);

    gr_ctx_init_nmod8(ZZn, 10);
    gr_test_ring(ZZn, 10000, flags);
    gr_ctx_clear(ZZn);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
