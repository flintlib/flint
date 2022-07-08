#include "gr.h"

int main()
{
    gr_ctx_t ZZ;
    int flags = 0;

    flint_printf("fmpz....");
    fflush(stdout);

    gr_ctx_init_fmpz(ZZ);
    ZZ->size_limit = 1000;
    gr_test_ring(ZZ, 10000, flags);
    gr_ctx_clear(ZZ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
