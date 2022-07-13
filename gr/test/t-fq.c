#include "gr.h"

int main()
{
    gr_ctx_t Fq;
    fmpz_t p;
    slong d;
    slong iter;
    flint_rand_t state;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fq....");
    fflush(stdout);

    flint_randinit(state);
    fmpz_init(p);

    for (iter = 0; iter < 30; iter++)
    {
        fmpz_randprime(p, state, 2 + n_randint(state, 100), 0);
        d = 1 + n_randint(state, 5);
        gr_ctx_init_fq(Fq, p, d, "a");
        gr_test_ring(Fq, 100, flags);
        gr_ctx_clear(Fq);
    }

    fmpz_clear(p);
    flint_randclear(state);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
