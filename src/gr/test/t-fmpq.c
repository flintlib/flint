#include "gr.h"

int main(void)
{
    gr_ctx_t QQ;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fmpq....");
    fflush(stdout);

    gr_ctx_init_fmpq(QQ);
    QQ->size_limit = 1000;
    gr_test_ring(QQ, 10000, flags);
    gr_ctx_clear(QQ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
