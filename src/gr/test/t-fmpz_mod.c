#include "fmpz.h"
#include "gr.h"

int main(void)
{
    gr_ctx_t ZZn;
    fmpz_t n;
    slong iter;
    flint_rand_t state;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fmpz_mod....");
    fflush(stdout);

    flint_randinit(state);
    fmpz_init(n);

    for (iter = 0; iter < 100; iter++)
    {
        fmpz_randtest_not_zero(n, state, 200);
        fmpz_abs(n, n);
        gr_ctx_init_fmpz_mod(ZZn, n);
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    /* test huge preinvn code */
    {
        fmpz_randbits(n, state, 72000);
        fmpz_abs(n, n);
        gr_ctx_init_fmpz_mod(ZZn, n);
        gr_test_ring(ZZn, 10, flags);
        gr_ctx_clear(ZZn);
    }

    fmpz_clear(n);
    flint_randclear(state);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
