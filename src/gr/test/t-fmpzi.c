#include "gr.h"

int main(void)
{
    gr_ctx_t ZZi;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fmpzi....");
    fflush(stdout);

    gr_ctx_init_fmpzi(ZZi);
    ZZi->size_limit = 1000;
    gr_test_ring(ZZi, 100000, flags);

    gr_ctx_clear(ZZi);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
