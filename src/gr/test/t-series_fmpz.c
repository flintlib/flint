#include "gr.h"

int main(void)
{
    gr_ctx_t ZZ, ZZx;
    int flags = 0;
    slong i;

    flint_printf("series_fmpz....");
    fflush(stdout);

    gr_ctx_init_fmpz(ZZ);

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_gr_series(ZZx, ZZ, i);
        gr_test_ring(ZZx, 100, flags);
        gr_ctx_clear(ZZx);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_gr_series_mod(ZZx, ZZ, i);
        gr_test_ring(ZZx, 100, flags);
        gr_ctx_clear(ZZx);
    }

    gr_ctx_clear(ZZ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
