
#include "gr_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("invert_rows_cols...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100; iter++)
    {
        slong m, n, i, j;
        gr_ctx_t ctx;
        gr_mat_t A;
        gr_mat_t B;

        gr_ctx_init_random(ctx, state);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(B, m, n, ctx);

        if (gr_mat_randtest(A, state, ctx) == GR_SUCCESS &&
            gr_mat_set(B, A, ctx) == GR_SUCCESS)
        {
            GR_MUST_SUCCEED(gr_mat_invert_rows(A, NULL, ctx));
            GR_MUST_SUCCEED(gr_mat_invert_cols(A, NULL, ctx));

            for (i = 0; i < A->r; i++)
            {
                for (j = 0; j < A->c; j++)
                {
                    if (gr_equal(gr_mat_entry_ptr(B, i, j, ctx), gr_mat_entry_ptr(A, A->r - i - 1, A->c - j - 1, ctx), ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: B != A\n");
                        flint_printf("A:\n");
                        gr_mat_print(A, ctx);
                        flint_printf("B:\n");
                        gr_mat_print(B, ctx);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
