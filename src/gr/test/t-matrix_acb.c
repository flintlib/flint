#include "gr.h"

int main(void)
{
    gr_ctx_t CC, MatCC;
    slong n;
    int flags = 0;

    flint_printf("matrix_acb....");
    fflush(stdout);

    gr_ctx_init_complex_acb(CC, 64);

    for (n = 0; n <= 5; n++)
    {
        gr_ctx_init_matrix_ring(MatCC, CC, n);
        gr_test_ring(MatCC, n <= 2 ? 100 : 10, flags);
        gr_ctx_clear(MatCC);
    }

    gr_ctx_clear(CC);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
