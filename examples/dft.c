/* This file is public domain. Author: Fredrik Johansson. */

#include <stdlib.h>
#include <string.h>
#include <flint/profiler.h>
#include <flint/calcium.h>
#include <flint/ca.h>
#include <flint/ca_vec.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_special.h>
#include <flint/gr_tower.h>
#include <flint/gr_tower_lazy.h>

void
benchmark_DFT(slong N, int input, int verbose, slong qqbar_limit, slong gb, ca_ctx_t ctx)
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
    /* ctx->options[CA_OPT_VERBOSE] = 1; */

    ctx->options[CA_OPT_USE_GROEBNER] = gb;

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
        else if (input == 6 || input == 7)
        {
            /* (from an fmpz: ca_pow_ui would create a symbolic power
               beyond CA_OPT_POW_LIMIT) */
            fmpz_t c;
            fmpz_init(c);
            fmpz_ui_pow_ui(c, i + 2, (input == 6) ? 1000 : 10000);
            ca_set_fmpz(x + i, c, ctx);
            fmpz_clear(c);
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

/* The same benchmark for a generic gr context (used with the lazy tower
   field, gr_ctx_init_tower_lazy). */
void
benchmark_DFT_gr(slong N, int input, int verbose, gr_ctx_t ctx)
{
    gr_ptr x, X, y, w, t;
    slong i, k, n, sz = ctx->sizeof_elem;
    truth_t is_zero;

#define E(v, i) GR_ENTRY(v, i, sz)

    x = gr_heap_init_vec(N, ctx);
    X = gr_heap_init_vec(N, ctx);
    y = gr_heap_init_vec(N, ctx);
    w = gr_heap_init_vec(2 * N, ctx);
    t = gr_heap_init(ctx);

    if (verbose)
        flint_printf("[x] =\n");

    for (i = 0; i < N; i++)
    {
        if (input == 0)
        {
            GR_MUST_SUCCEED(gr_set_ui(E(x, i), i + 2, ctx));
        }
        else if (input == 1)
        {
            GR_MUST_SUCCEED(gr_set_ui(E(x, i), i + 2, ctx));
            GR_MUST_SUCCEED(gr_sqrt(E(x, i), E(x, i), ctx));
        }
        else if (input == 2)
        {
            GR_MUST_SUCCEED(gr_set_ui(E(x, i), i + 2, ctx));
            GR_MUST_SUCCEED(gr_log(E(x, i), E(x, i), ctx));
        }
        else if (input == 3)
        {
            GR_MUST_SUCCEED(gr_pi(E(x, i), ctx));
            GR_MUST_SUCCEED(gr_i(t, ctx));
            GR_MUST_SUCCEED(gr_mul(E(x, i), E(x, i), t, ctx));
            GR_MUST_SUCCEED(gr_mul_ui(E(x, i), E(x, i), 2, ctx));
            GR_MUST_SUCCEED(gr_div_ui(E(x, i), E(x, i), i + 2, ctx));
            GR_MUST_SUCCEED(gr_exp(E(x, i), E(x, i), ctx));
        }
        else if (input == 4)
        {
            GR_MUST_SUCCEED(gr_pi(E(x, i), ctx));
            GR_MUST_SUCCEED(gr_mul_ui(E(x, i), E(x, i), i + 2, ctx));
            GR_MUST_SUCCEED(gr_add_ui(E(x, i), E(x, i), 1, ctx));
            GR_MUST_SUCCEED(gr_inv(E(x, i), E(x, i), ctx));
        }
        else if (input == 5)
        {
            GR_MUST_SUCCEED(gr_pi(E(x, i), ctx));
            GR_MUST_SUCCEED(gr_set_ui(t, i + 2, ctx));
            GR_MUST_SUCCEED(gr_sqrt(t, t, ctx));
            GR_MUST_SUCCEED(gr_mul(E(x, i), E(x, i), t, ctx));
            GR_MUST_SUCCEED(gr_add_ui(E(x, i), E(x, i), 1, ctx));
            GR_MUST_SUCCEED(gr_inv(E(x, i), E(x, i), ctx));
        }
        else if (input == 6 || input == 7)
        {
            fmpz_t c;
            fmpz_init(c);
            fmpz_ui_pow_ui(c, i + 2, (input == 6) ? 1000 : 10000);
            GR_MUST_SUCCEED(gr_set_fmpz(E(x, i), c, ctx));
            fmpz_clear(c);
        }

        if (verbose)
            gr_println(E(x, i), ctx);
    }

    /* roots of unity */
    for (i = 0; i < 2 * N; i++)
    {
        if (i == 0)
        {
            GR_MUST_SUCCEED(gr_one(E(w, i), ctx));
        }
        else if (i == 1)
        {
            GR_MUST_SUCCEED(gr_pi(E(w, i), ctx));
            GR_MUST_SUCCEED(gr_i(t, ctx));
            GR_MUST_SUCCEED(gr_mul(E(w, i), E(w, i), t, ctx));
            GR_MUST_SUCCEED(gr_mul_ui(E(w, i), E(w, i), 2, ctx));
            GR_MUST_SUCCEED(gr_div_si(E(w, i), E(w, i), N, ctx));
            GR_MUST_SUCCEED(gr_exp(E(w, i), E(w, i), ctx));
        }
        else
        {
            GR_MUST_SUCCEED(gr_mul(E(w, i), E(w, i - 1), E(w, 1), ctx));
        }
    }

    if (verbose)
        printf("\nDFT([x]) =\n");

    for (k = 0; k < N; k++)
    {
        GR_MUST_SUCCEED(gr_zero(E(X, k), ctx));
        for (n = 0; n < N; n++)
        {
            GR_MUST_SUCCEED(gr_mul(t, E(x, n), E(w, ((2 * N - k) * n) % (2 * N)), ctx));
            GR_MUST_SUCCEED(gr_add(E(X, k), E(X, k), t, ctx));
        }
        if (verbose)
            gr_println(E(X, k), ctx);
    }

    if (verbose)
        printf("\nIDFT(DFT([x])) =\n");

    for (k = 0; k < N; k++)
    {
        GR_MUST_SUCCEED(gr_zero(E(y, k), ctx));
        for (n = 0; n < N; n++)
        {
            GR_MUST_SUCCEED(gr_mul(t, E(X, n), E(w, (k * n) % (2 * N)), ctx));
            GR_MUST_SUCCEED(gr_add(E(y, k), E(y, k), t, ctx));
        }
        GR_MUST_SUCCEED(gr_div_ui(E(y, k), E(y, k), N, ctx));
        if (verbose)
            gr_println(E(y, k), ctx);
    }

    if (verbose)
        printf("\n[x] - IDFT(DFT([x])) =\n");

    for (k = 0; k < N; k++)
    {
        GR_MUST_SUCCEED(gr_sub(t, E(x, k), E(y, k), ctx));
        is_zero = gr_is_zero(t, ctx);
        if (verbose)
        {
            gr_print(t, ctx);
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

    gr_heap_clear_vec(x, N, ctx);
    gr_heap_clear_vec(X, N, ctx);
    gr_heap_clear_vec(y, N, ctx);
    gr_heap_clear_vec(w, 2 * N, ctx);
    gr_heap_clear(t, ctx);
#undef E
}

void usage(void)
{
    printf("usage: dft [-verbose] [-input i] [-limit B] [-timing T] [-nogb] [-tower] [-gens F] [-cyclo D] N\n");
}

int main(int argc, char *argv[])
{
    ca_ctx_t ctx;
    int verbose, input, timing;
    slong i, Nmin, Nmax, N, qqbar_limit, gb;
    int tower = 0, gens = 0;
    slong cyclo = -1;

    Nmin = Nmax = 2;
    verbose = 0;
    input = 0;
    timing = 0;
    qqbar_limit = 0;
    gb = 1;

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
        else if (!strcmp(argv[i], "-nogb"))
        {
            gb = 0;
        }
        else if (!strcmp(argv[i], "-tower"))
        {
            tower = 1;
        }
        else if (!strcmp(argv[i], "-cyclo"))
        {
            /* the option GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT (implies -tower) */
            tower = 1;
            cyclo = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-gens"))
        {
            /* generator policy flags of the lazy tower field
               (GR_TOWER_GENS_*: 2 = composite roots of unity, 4 =
               composite square roots) */
            gens = atol(argv[i+1]);
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
        flint_printf("DFT benchmark, length N = %wd\n", N);
        if (input == 0)
            flint_printf("x_k = k + 2\n");
        else if (input == 1)
            flint_printf("x_k = sqrt(k + 2)\n");
        else if (input == 2)
            flint_printf("x_k = log(k + 2)\n");
        else if (input == 3)
            flint_printf("x_k = exp(2 pi i / (k + 2))\n");
        else if (input == 4)
            flint_printf("x_k = 1 / (1 + (k + 2) pi)\n");
        else if (input == 5)
            flint_printf("x_k = 1 / (1 + sqrt(k + 2) pi)\n");
        else if (input == 6)
            flint_printf("x_k = (k + 2)^1000\n");
        else if (input == 7)
            flint_printf("x_k = (k + 2)^10000\n");

        flint_printf("\n");

        if (tower)
        {
            gr_ctx_t QQ, K;
            TIMEIT_ONCE_START;
            gr_ctx_init_fmpq(QQ);
            gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
            gr_tower_lazy_ctx_set_gen_flags(K, gens);
            if (cyclo >= 0)
                gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT, cyclo);
            benchmark_DFT_gr(N, input, verbose, K);
            gr_tower_lazy_ctx_stats(K);
            gr_ctx_clear(K);
            gr_ctx_clear(QQ);
            TIMEIT_ONCE_STOP;
        }
        else if (timing == 0)
        {
            TIMEIT_ONCE_START;
            ca_ctx_init(ctx);
            benchmark_DFT(N, input, verbose, qqbar_limit, gb, ctx);
            ca_ctx_clear(ctx);
            TIMEIT_ONCE_STOP;
        }
        else if (timing == 1)
        {
            TIMEIT_START;
            ca_ctx_init(ctx);
            benchmark_DFT(N, input, verbose, qqbar_limit, gb, ctx);
            ca_ctx_clear(ctx);
            TIMEIT_STOP;
        }
        else
        {
            ca_ctx_init(ctx);
            benchmark_DFT(N, input, verbose, qqbar_limit, gb, ctx);
            TIMEIT_START;
            benchmark_DFT(N, input, verbose, qqbar_limit, gb, ctx);
            TIMEIT_STOP;
            ca_ctx_clear(ctx);
        }
    }

    print_memory_usage();
    flint_cleanup();
    return 0;
}
