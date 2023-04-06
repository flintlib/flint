#include "gr.h"

int main(void)
{
    gr_ctx_t RR, RRx;
    int flags = 0;

    flint_printf("polynomial_arb....");
    fflush(stdout);

    gr_ctx_init_real_arb(RR, 64);
    gr_ctx_init_polynomial(RRx, RR);
    RRx->size_limit = 50;
    gr_test_ring(RRx, 1000, flags);
    gr_ctx_clear(RRx);
    gr_ctx_clear(RR);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
