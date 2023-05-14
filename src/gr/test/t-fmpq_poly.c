#include "gr.h"

int main(void)
{
    gr_ctx_t QQx;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fmpq_poly....");
    fflush(stdout);

    gr_ctx_init_fmpq_poly(QQx);
    gr_test_ring(QQx, 1000, flags);
    gr_ctx_clear(QQx);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
