#include "gr.h"

int main(void)
{
    gr_ctx_t RR, CC, QQbar_real, QQbar;
    int flags = 0;

    flint_printf("ca....");
    fflush(stdout);

    gr_ctx_init_real_ca(RR);
    gr_test_ring(RR, 100, flags);
    gr_ctx_clear(RR);

    gr_ctx_init_complex_ca(CC);
    gr_test_ring(CC, 100, flags);
    gr_ctx_clear(CC);

    gr_ctx_init_real_algebraic_ca(QQbar_real);
    gr_test_ring(QQbar_real, 100, flags);
    gr_ctx_clear(QQbar_real);

    gr_ctx_init_complex_algebraic_ca(QQbar);
    gr_test_ring(QQbar, 100, flags);
    gr_ctx_clear(QQbar);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
