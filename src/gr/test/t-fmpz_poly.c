#include "gr.h"

int main()
{
    gr_ctx_t ZZx;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fmpz_poly....");
    fflush(stdout);

    gr_ctx_init_fmpz_poly(ZZx);
    gr_test_ring(ZZx, 1000, flags);

    gr_ctx_clear(ZZx);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
