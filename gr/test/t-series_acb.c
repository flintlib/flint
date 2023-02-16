#include "gr.h"

int main()
{
    gr_ctx_t CCn, CCnx;
    int flags = 0;
    slong i;
    flint_rand_t state;

    flint_randinit(state);

    flint_printf("series_acb....");
    fflush(stdout);

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_complex_acb(CCn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series(CCnx, CCn, i);
        gr_test_ring(CCnx, 100, flags);
        gr_ctx_clear(CCnx);
        gr_ctx_clear(CCn);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_complex_acb(CCn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series_mod(CCnx, CCn, i);
        gr_test_ring(CCnx, 100, flags);
        gr_ctx_clear(CCnx);
        gr_ctx_clear(CCn);
    }

    gr_ctx_clear(CCn);

    flint_randclear(state);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
