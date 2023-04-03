
#include "gr_mat.h"

/* todo: make a function */
int
gr_mat_randrowops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (mat->r == 0 || mat->c == 0)
        return GR_SUCCESS;

    for (c = 0; c < count; c++)
    {
        if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
            continue;
        if (n_randint(state, 2))
            for (k = 0; k < n; k++)
                status |= gr_add(GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, i, k, sz), ctx);
        else
            for (k = 0; k < n; k++)
                status |= gr_sub(GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, i, k, sz), ctx);
    }

    return status;
}

int main()
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;
    flint_rand_t state;

    flint_printf("rref_fflu...");
    fflush(stdout);

    flint_randinit(state);

    /* Check that random row/column operations preserve rank */
    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, R, R2;
        slong rank1, rank2, r, c;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(B, r, c, ctx);
        gr_mat_init(R, r, c, ctx);
        gr_mat_init(R2, r, c, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_set(B, A, ctx);
        status |= gr_mat_randrowops(B, state, 1 + n_randint(state, 20), ctx);

        status |= gr_mat_rref_fflu(&rank1, R, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_mat_rref_fflu(&rank2, R2, B, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank1 != rank2 || gr_mat_equal(R, R2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("R: "); gr_mat_print(R, ctx); flint_printf("\n");
                    flint_printf("R2: "); gr_mat_print(R2, ctx); flint_printf("\n");
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
        gr_mat_clear(R, ctx);
        gr_mat_clear(R2, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" [%wd success, %wd domain, %wd unable] PASS\n", count_success, count_domain, count_unable);
    return 0;
}
