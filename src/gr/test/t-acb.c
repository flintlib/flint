#include "gr.h"

int main(void)
{
    gr_ctx_t CC;
    int flags = 0;
    slong prec;

    flint_printf("acb....");
    fflush(stdout);

    for (prec = 64; prec <= 256; prec *= 2)
    {
        gr_ctx_init_complex_acb(CC, prec);
        gr_test_ring(CC, 1000, flags);
        gr_ctx_clear(CC);
    }

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
