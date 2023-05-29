#include "mpoly.h"
#include "gr.h"

int main(void)
{
    gr_ctx_t ZZxy;
    slong iter;
    int flags = 0;
    flint_rand_t state;

    flint_printf("fmpz_mpoly....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10; iter++)
    {
        gr_ctx_init_fmpz_mpoly(ZZxy, n_randint(state, 3), mpoly_ordering_randtest(state));
        ZZxy->size_limit = 100;
        gr_test_ring(ZZxy, 100, flags);
        gr_ctx_clear(ZZxy);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
