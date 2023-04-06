
#include "gr_mat.h"

int main(void)
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;
    flint_rand_t state;

    flint_printf("hadamard...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        slong n;
        int status = GR_SUCCESS;

        if (n_randint(state, 20) == 0)
        {
            gr_ctx_init_fmpz(ctx);
            n = 4 * n_randint(state, 30);
        }
        else
        {
            gr_ctx_init_random(ctx, state);
            n = n_randint(state, 20);
        }

        gr_mat_init(A, n, n, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        status = gr_mat_hadamard(A, ctx);

        if (status == GR_SUCCESS)
        {
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            status |= gr_mat_transpose(B, A, ctx);
            status |= gr_mat_mul(C, A, B, ctx);
            status |= gr_mat_set_ui(B, n, ctx);

            if (status == GR_SUCCESS && gr_mat_equal(C, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                gr_ctx_println(ctx);
                flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                flint_printf("C: "); gr_mat_print(C, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" [%wd success, %wd domain, %wd unable] PASS\n", count_success, count_domain, count_unable);
    return 0;
}
