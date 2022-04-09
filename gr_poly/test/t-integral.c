#include "gr_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("integral....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t F, G, Gprime;

        /* check derivative(integral(f)) == f */
        gr_ctx_init_random(ctx, state);

        gr_poly_init(F, ctx);
        gr_poly_init(G, ctx);
        gr_poly_init(Gprime, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(F, state, 1 + n_randint(state, 6), ctx);
        status |= gr_poly_randtest(G, state, 1 + n_randint(state, 6), ctx);
        status |= gr_poly_randtest(Gprime, state, 1 + n_randint(state, 6), ctx);

        if (n_randint(state, 2))
        {
            status |= gr_poly_integral(G, F, ctx);
        }
        else
        {
            status |= gr_poly_set(G, F, ctx);
            status |= gr_poly_integral(G, G, ctx);
        }

        if (n_randint(state, 2))
        {
            status |= gr_poly_derivative(Gprime, G, ctx);
        }
        else
        {
            status |= gr_poly_set(Gprime, G, ctx);
            status |= gr_poly_derivative(Gprime, Gprime, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_equal(F, Gprime, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
            flint_printf("Gprime = "); gr_poly_print(Gprime, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(F, ctx);
        gr_poly_clear(G, ctx);
        gr_poly_clear(Gprime, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
