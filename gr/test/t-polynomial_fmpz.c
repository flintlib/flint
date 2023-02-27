#include "gr.h"

int main()
{
    gr_ctx_t ZZ, ZZx;
    int flags = 0;

    flint_printf("polynomial_fmpz....");
    fflush(stdout);

    gr_ctx_init_fmpz(ZZ);
    ZZ->size_limit = 200;
    gr_ctx_init_polynomial(ZZx, ZZ);
    ZZx->size_limit = 30;
    gr_test_ring(ZZx, 1000, flags);
    gr_ctx_clear(ZZx);
    gr_ctx_clear(ZZ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
