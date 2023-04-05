
#include "gr_mat.h"

int main(void)
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;
    flint_rand_t state;

    flint_printf("rank...");
    fflush(stdout);

    flint_randinit(state);

    /* Check that random row/column operations preserve rank */
    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B;
        slong rank1, rank2, r, c;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(B, r, c, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_set(B, A, ctx);
        status |= gr_mat_randops(B, state, 1 + n_randint(state, 20), ctx);

        status |= gr_mat_rank(&rank1, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_mat_rank(&rank2, B, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank1 != rank2)
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" [%wd success, %wd domain, %wd unable] PASS\n", count_success, count_domain, count_unable);
    return 0;
}
