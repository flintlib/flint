/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "ca.h"

int main(int argc, char *argv[])
{
    ca_ctx_t ctx;
    ca_t sqrt5, phi, psi, t, u;
    fmpz_t n;

    if (argc < 2)
    {
        flint_printf("usage: build/examples/binet [-limit B] n\n");
        return 1;
    }

    fmpz_init(n);
    fmpz_set_str(n, argv[argc-1], 10);

    TIMEIT_ONCE_START
    ca_ctx_init(ctx);

    if (argc == 4)
        ctx->options[CA_OPT_PREC_LIMIT] = atol(argv[2]);

    ca_init(sqrt5, ctx);
    ca_init(phi, ctx);
    ca_init(psi, ctx);
    ca_init(t, ctx);
    ca_init(u, ctx);

    ca_sqrt_ui(sqrt5, 5, ctx);

    ca_add_ui(phi, sqrt5, 1, ctx);
    ca_div_ui(phi, phi, 2, ctx);

    ca_ui_sub(psi, 1, phi, ctx);

    ca_pow_fmpz(t, phi, n, ctx);
    ca_pow_fmpz(u, psi, n, ctx);

    ca_sub(t, t, u, ctx);
    ca_div(t, t, sqrt5, ctx);

    ca_print(t, ctx);
    flint_printf("\n");

    ca_clear(sqrt5, ctx);
    ca_clear(phi, ctx);
    ca_clear(psi, ctx);
    ca_clear(t, ctx);
    ca_clear(u, ctx);
    ca_ctx_clear(ctx);

    fmpz_clear(n);
    flint_printf("\n");
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}
