
#include "gr_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("randrank...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        slong rank1, rank2, r, c;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        if (gr_ctx_is_integral_domain(ctx) == T_TRUE)
        {
            r = n_randint(state, 6);
            c = n_randint(state, 6);

            gr_mat_init(A, r, c, ctx);

            rank1 = n_randint(state, 6);

            status |= gr_mat_randrank(A, state, rank1, ctx);
            status |= gr_mat_rank(&rank2, A, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank1 != rank2)
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }

            gr_mat_clear(A, ctx);
        }

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" PASS\n");
    return EXIT_SUCCESS;
}
