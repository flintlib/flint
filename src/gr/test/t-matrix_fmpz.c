#include "gr.h"

int main(void)
{
    gr_ctx_t ZZ, MatZZ;
    slong n;
    int flags = 0;

    flint_printf("matrix_fmpz....");
    fflush(stdout);

    gr_ctx_init_fmpz(ZZ);
    ZZ->size_limit = 1000;

    for (n = 0; n <= 8; n++)
    {
        gr_ctx_init_matrix_ring(MatZZ, ZZ, n);
        gr_test_ring(MatZZ, n <= 2 ? 100 : 10, flags);
        gr_ctx_clear(MatZZ);
    }

    gr_ctx_clear(ZZ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
