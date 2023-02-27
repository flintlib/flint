#include "gr.h"

int main()
{
    gr_ctx_t RR, MatRR;
    slong n;
    int flags = 0;

    flint_printf("matrix_arb....");
    fflush(stdout);

    gr_ctx_init_real_arb(RR, 64);

    for (n = 0; n <= 5; n++)
    {
        gr_ctx_init_matrix_ring(MatRR, RR, n);
        gr_test_ring(MatRR, n <= 2 ? 100 : 10, flags);
        gr_ctx_clear(MatRR);
    }

    gr_ctx_clear(RR);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
