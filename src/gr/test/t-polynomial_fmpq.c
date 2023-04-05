#include "gr.h"

int main(void)
{
    gr_ctx_t QQ, QQx;
    int flags = 0;

    flint_printf("polynomial_fmpq....");
    fflush(stdout);

    gr_ctx_init_fmpq(QQ);
    QQ->size_limit = 100;
    gr_ctx_init_polynomial(QQx, QQ);
    QQx->size_limit = 30;
    gr_test_ring(QQx, 1000, flags);
    gr_ctx_clear(QQx);
    gr_ctx_clear(QQ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
