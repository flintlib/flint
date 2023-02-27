#include "gr.h"

int main()
{
    gr_ctx_t RR;
    int flags = 0;
    slong prec;

    flint_printf("arb....");
    fflush(stdout);

    for (prec = 64; prec <= 256; prec *= 2)
    {
        gr_ctx_init_real_arb(RR, prec);
        gr_test_ring(RR, 1000, flags);
        gr_ctx_clear(RR);
    }

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
