#include "gr.h"

int main(void)
{
    gr_ctx_t QQ, MatQQ;
    slong n;
    int flags = 0;

    flint_printf("matrix_fmpq....");
    fflush(stdout);

    gr_ctx_init_fmpq(QQ);
    QQ->size_limit = 200;

    for (n = 0; n <= 5; n++)
    {
        gr_ctx_init_matrix_ring(MatQQ, QQ, n);
        gr_test_ring(MatQQ, n <= 2 ? 100 : 10, flags);
        gr_ctx_clear(MatQQ);
    }

    gr_ctx_clear(QQ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
