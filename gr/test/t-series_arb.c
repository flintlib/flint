#include "flint/ulong_extras.h"
#include "gr.h"

int main()
{
    gr_ctx_t RRn, RRnx;
    int flags = 0;
    slong i;
    flint_rand_t state;

    flint_randinit(state);

    flint_printf("series_arb....");
    fflush(stdout);

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_real_arb(RRn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series(RRnx, RRn, i);
        gr_test_ring(RRnx, 100, flags);
        gr_ctx_clear(RRnx);
        gr_ctx_clear(RRn);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_real_arb(RRn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series_mod(RRnx, RRn, i);
        gr_test_ring(RRnx, 100, flags);
        gr_ctx_clear(RRnx);
        gr_ctx_clear(RRn);
    }

    gr_ctx_clear(RRn);

    flint_randclear(state);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
