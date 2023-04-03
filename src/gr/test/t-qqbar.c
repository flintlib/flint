#include "gr.h"

int main()
{
    gr_ctx_t QQbar_real, QQbar;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("qqbar....");
    fflush(stdout);

    gr_ctx_init_real_qqbar(QQbar_real);
    gr_test_ring(QQbar_real, 100, flags);
    gr_ctx_clear(QQbar_real);

    gr_ctx_init_complex_qqbar(QQbar);
    gr_test_ring(QQbar, 100, flags);
    gr_ctx_clear(QQbar);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
