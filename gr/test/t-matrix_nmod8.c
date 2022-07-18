#include "gr.h"

int main()
{
    gr_ctx_t ZZn, MatZZn;
    slong iter, n;
    int flags = 0;
    flint_rand_t state;

    flint_printf("matrix_nmod8....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 50; iter++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);

        for (n = 0; n <= 4; n++)
        {
            gr_ctx_init_matrix_ring(MatZZn, ZZn, n);
            gr_test_ring(MatZZn, n <= 2 ? 100 : 10, flags);
            gr_ctx_clear(MatZZn);
        }

        if (iter % 5 == 0)
        {
            gr_ctx_init_matrix_ring(MatZZn, ZZn, n_randint(state, 40));
            gr_test_ring(MatZZn, 10, flags);
            gr_ctx_clear(MatZZn);
        }

        gr_ctx_clear(ZZn);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
