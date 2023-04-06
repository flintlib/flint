#include "gr.h"

int main(void)
{
    gr_ctx_t CC, CCx;
    int flags = 0;

    flint_printf("polynomial_acb....");
    fflush(stdout);

    gr_ctx_init_complex_acb(CC, 64);
    gr_ctx_init_polynomial(CCx, CC);
    CCx->size_limit = 50;
    gr_test_ring(CCx, 1000, flags);
    gr_ctx_clear(CCx);
    gr_ctx_clear(CC);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
