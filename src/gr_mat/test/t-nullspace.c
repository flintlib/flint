
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
        gr_mat_t A, X, AX;
        slong r, c, rank, nullity;
        int status = GR_SUCCESS;
        int status2 = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(X, 0, 0, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        status = gr_mat_nullspace(X, A, ctx);

        if (status == GR_SUCCESS)
        {
            status2 = gr_mat_rank(&rank, A, ctx);

            if (status2 == GR_SUCCESS)
            {
                nullity = c - rank;

                if (nullity != gr_mat_ncols(X, ctx))
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("r = %wd, c = %wd\n", r, c);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); gr_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, nullity = %wd\n\n", rank, nullity);
                    flint_abort();
                }

                gr_mat_init(AX, r, nullity, ctx);
                status2 |= gr_mat_mul(AX, A, X, ctx);

                if (status2 == GR_SUCCESS && gr_mat_is_zero(AX, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (2):\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); gr_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("AX: "); gr_mat_print(AX, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, nullity = %wd\n\n", rank, nullity);
                    flint_abort();
                }

                status2 |= gr_mat_rank(&nullity, X, ctx);

                if (status2 == GR_SUCCESS && nullity != gr_mat_ncols(X, ctx))
                {
                    flint_printf("FAIL (3):\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); gr_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank, nullity);
                    flint_abort();
                }

                gr_mat_clear(AX, ctx);
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(X, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" [%wd success, %wd domain, %wd unable] PASS\n", count_success, count_domain, count_unable);
    return 0;
}
