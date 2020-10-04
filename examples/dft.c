/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "ca.h"
#include "ca_vec.h"

void
benchmark_DFT(slong N, int input, int verbose, slong qqbar_limit, ca_ctx_t ctx)
{
    ca_ptr x, X, y, w;
    ca_t t;
    slong i, k, n;
    truth_t is_zero;

    x = _ca_vec_init(N, ctx);
    X = _ca_vec_init(N, ctx);
    y = _ca_vec_init(N, ctx);
    w = _ca_vec_init(2 * N, ctx);
    ca_init(t, ctx);

    /* ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_DEBUG; */

    if (qqbar_limit != 0)
        ctx->options[CA_OPT_QQBAR_DEG_LIMIT] = qqbar_limit;

    /* Construct input vector */
    if (verbose)
        flint_printf("[x] =\n");
    for (i = 0; i < N; i++)
    {
        if (input == 0)
        {
            ca_set_ui(x + i, i + 2, ctx);
        }
        else if (input == 1)
        {
            ca_set_ui(x + i, i + 2, ctx);
            ca_sqrt(x + i, x + i, ctx);
        }
        else if (input == 2)
        {
            ca_set_ui(x + i, i + 2, ctx);
            ca_log(x + i, x + i, ctx);
        }
        else if (input == 3)
        {
            ca_pi_i(x + i, ctx);
            ca_mul_ui(x + i, x + i, 2, ctx);
            ca_div_ui(x + i, x + i, i + 2, ctx);
            ca_exp(x + i, x + i, ctx);
        }
        else if (input == 4)
        {
            ca_pi(x + i, ctx);
            ca_mul_ui(x + i, x + i, i + 2, ctx);
            ca_add_ui(x + i, x + i, 1, ctx);
            ca_inv(x + i, x + i, ctx);
        }
        else if (input == 5)
        {
            ca_pi(x + i, ctx);
            ca_sqrt_ui(w, i + 2, ctx);
            ca_mul(x + i, x + i, w, ctx);
            ca_add_ui(x + i, x + i, 1, ctx);
            ca_inv(x + i, x + i, ctx);
        }

        if (verbose)
        {
            ca_print(x + i, ctx);
            printf("\n");
        }
    }

    /* Construct roots of unity */
    for (i = 0; i < 2 * N; i++)
    {
        if (i == 0)
        {
            ca_one(w + i, ctx);
        }
        else if (i == 1)
        {
            ca_pi_i(w + i, ctx);
            ca_mul_ui(w + i, w + i, 2, ctx);
            ca_div_si(w + i, w + i, N, ctx);
            ca_exp(w + i, w + i, ctx);
        }
        else
        {
            ca_mul(w + i, w + i - 1, w + 1, ctx);
        }
    }

    /* Forward DFT */
    if (verbose)
        printf("\nDFT([x]) =\n");
    for (k = 0; k < N; k++)
    {
        ca_zero(X + k, ctx);

        for (n = 0; n < N; n++)
        {
            ca_mul(t, x + n, w + ((2 * N - k) * n) % (2 * N), ctx);
            ca_add(X + k, X + k, t, ctx);
        }

        if (verbose)
        {
            ca_print(X + k, ctx);
            printf("\n");
        }
    }

    /* Inverse DFT */
    if (verbose)
        printf("\nIDFT(DFT([x])) =\n");
    for (k = 0; k < N; k++)
    {
        ca_zero(y + k, ctx);

        for (n = 0; n < N; n++)
        {
            ca_mul(t, X + n, w + (k * n) % (2 * N), ctx);
            ca_add(y + k, y + k, t, ctx);

        }

        ca_div_ui(y + k, y + k, N, ctx);

        if (verbose)
        {
            ca_print(y + k, ctx);
            flint_printf("\n");
        }
    }

    if (verbose)
        printf("\n[x] - IDFT(DFT([x])) =\n");
    for (k = 0; k < N; k++)
    {
        ca_sub(t, x + k, y + k, ctx);

        is_zero = ca_check_is_zero(t, ctx);

        if (verbose)
        {
            ca_print(t, ctx);
            printf("       (= 0   ");
            truth_print(is_zero);
            printf(")\n");
        }

        if (is_zero != T_TRUE)
        {
            printf("Failed to prove equality!\n");
            flint_abort();
        }
    }

    if (verbose)
        printf("\n");

    _ca_vec_clear(x, N, ctx);
    _ca_vec_clear(X, N, ctx);
    _ca_vec_clear(y, N, ctx);
    _ca_vec_clear(w, 2 * N, ctx);
    ca_clear(t, ctx);
}

void usage()
{
    printf("usage: dft [-verbose] [-input i] [-limit B] [-timing T] N\n");
}

int main(int argc, char *argv[])
{
    ca_ctx_t ctx;
    int verbose, input, timing;
    slong i, Nmin, Nmax, N, qqbar_limit;

    Nmin = Nmax = 2;
    verbose = 0;
    input = 0;
    timing = 0;
    qqbar_limit = 0;

    if (argc < 2)
    {
        usage();
        return 1;
    }

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-verbose"))
        {
            verbose = 1;
        }
        else if (!strcmp(argv[i], "-input"))
        {
            input = atol(argv[i+1]);
            i += 1;
        }
        else if (!strcmp(argv[i], "-limit"))
        {
            qqbar_limit = atol(argv[i+1]);
            i += 1;
        }
        else if (!strcmp(argv[i], "-timing"))
        {
            timing = atol(argv[i+1]);
            i += 1;
        }
        else
        {
            Nmin = Nmax = atol(argv[i]);
            if (Nmin < 0)
            {
                Nmin = 0;
                Nmax = -Nmax;
            }
        }
    }

    for (N = Nmin; N <= Nmax; N++)
    {
        printf("DFT benchmark, length N = %ld\n\n", N);

        if (timing == 0)
        {
            TIMEIT_ONCE_START
            ca_ctx_init(ctx);
            benchmark_DFT(N, input, verbose, qqbar_limit, ctx);
            ca_ctx_clear(ctx);
            TIMEIT_ONCE_STOP
        }
        else if (timing == 1)
        {
            TIMEIT_START
            ca_ctx_init(ctx);
            benchmark_DFT(N, input, verbose, qqbar_limit, ctx);
            ca_ctx_clear(ctx);
            TIMEIT_STOP
        }
        else
        {
            ca_ctx_init(ctx);
            benchmark_DFT(N, input, verbose, qqbar_limit, ctx);
            TIMEIT_START
            benchmark_DFT(N, input, verbose, qqbar_limit, ctx);
            TIMEIT_STOP
            ca_ctx_clear(ctx);
        }
    }

    SHOW_MEMORY_USAGE
    flint_cleanup();
    return EXIT_SUCCESS;
}
