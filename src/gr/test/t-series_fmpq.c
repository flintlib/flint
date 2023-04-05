#include "gr.h"

int main(void)
{
    gr_ctx_t QQ, QQx;
    int flags = 0;
    slong i;

    flint_printf("series_fmpq....");
    fflush(stdout);

    gr_ctx_init_fmpq(QQ);

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_gr_series(QQx, QQ, i);
        gr_test_ring(QQx, 100, flags);
        gr_ctx_clear(QQx);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_gr_series_mod(QQx, QQ, i);
        gr_test_ring(QQx, 100, flags);
        gr_ctx_clear(QQx);
    }

    gr_ctx_clear(QQ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
