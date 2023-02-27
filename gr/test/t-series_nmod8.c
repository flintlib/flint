#include "ulong_extras.h"
#include "gr.h"

int main()
{
    gr_ctx_t ZZn, ZZnx;
    int flags = 0;
    slong i;
    flint_rand_t state;

    flint_randinit(state);

    flint_printf("series_nmod8....");
    fflush(stdout);

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);
        gr_ctx_init_gr_series(ZZnx, ZZn, i);
        gr_test_ring(ZZnx, 100, flags);
        gr_ctx_clear(ZZnx);
        gr_ctx_clear(ZZn);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);
        gr_ctx_init_gr_series_mod(ZZnx, ZZn, i);
        gr_test_ring(ZZnx, 100, flags);
        gr_ctx_clear(ZZnx);
        gr_ctx_clear(ZZn);
    }

    gr_ctx_clear(ZZn);

    flint_randclear(state);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
