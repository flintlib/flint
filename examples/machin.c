/* This file is public domain. Author: Fredrik Johansson. */

/* Verifies Machin-like formulas exactly with the ca (exact real and
   complex number) module: a few classical formulas for pi and for
   logarithms, and then the sets tabulated in FLINT itself, which are
   used for the precomputations behind the elementary functions.
   Usage: machin [maxn]   (default 8; the tabulated sets go to 48)
          machin -tower [maxn]   (the lazy tower field, gr_tower.h) */

#include <flint/profiler.h>
#include <stdlib.h>
#include <string.h>
#include <flint/ca.h>
#include <flint/ca_vec.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <flint/mp_real.h>
#include <flint/gr.h>
#include <flint/gr_special.h>
#include <flint/gr_tower.h>
#include <flint/gr_tower_lazy.h>

void
simple_ca_atan_p_q(ca_t res, ulong p, ulong q, ca_ctx_t ctx)
{
    ca_set_ui(res, p, ctx);
    ca_div_ui(res, res, q, ctx);
    ca_atan(res, res, ctx);
}

/* valid for -1 < x < 1 */
void
simple_ca_atanh(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t t, u;

    ca_init(t, ctx);
    ca_init(u, ctx);

    ca_add_ui(t, x, 1, ctx);
    ca_sub_ui(u, x, 1, ctx);
    ca_neg(u, u, ctx);
    ca_div(res, t, u, ctx);
    ca_log(res, res, ctx);
    ca_div_ui(res, res, 2, ctx);

    ca_clear(t, ctx);
    ca_clear(u, ctx);
}

void
simple_ca_atanh_p_q(ca_t res, ulong p, ulong q, ca_ctx_t ctx)
{
    ca_set_ui(res, p, ctx);
    ca_div_ui(res, res, q, ctx);
    simple_ca_atanh(res, res, ctx);
}

#define NUM_FORMULAS 8

slong machin_formulas[NUM_FORMULAS][4][2] = {
    {{1, 1}, {0, 0}, {0, 0}, {0, 0}},

    {{1, 2}, {1, 3}, {0, 0}, {0, 0}},
    {{2, 2}, {-1, 7}, {0, 0}, {0, 0}},
    {{2, 3}, {1, 7}, {0, 0}, {0, 0}},
    {{4, 5}, {-1, 239}, {0, 0}, {0, 0}},

    {{1, 2}, {1, 5}, {1, 8}, {0, 0}},
    {{1, 3}, {1, 4}, {1, 7}, {1, 13}},
    {{12, 49}, {32, 57}, {-5, 239}, {12, 110443}},
};

#define NUM_FORMULAS2 7

slong hyperbolic_logs[NUM_FORMULAS2] = {2, 3, 5, 2, 3, 5, 7};

slong hyperbolic_machin_formulas[NUM_FORMULAS2][4][2] = {
    {{14, 31}, {10, 49}, {6, 161}, {0, 0}},
    {{22, 31}, {16, 49}, {10, 161}, {0, 0}},
    {{32, 31}, {24, 49}, {14, 161}, {0, 0}},

    {{144, 251}, {54, 449}, {-38, 4801}, {62, 8749}},
    {{228, 251}, {86, 449}, {-60, 4801}, {98, 8749}},
    {{334, 251}, {126, 449}, {-88, 4801}, {144, 8749}},
    {{404, 251}, {152, 449}, {-106, 4801}, {174, 8749}},
};

/* Verify the Machin-type sets tabulated in FLINT (mp_real/machin_tab.c,
   _mp_real_machin_table): for the first n primes,

       den log(p_i) = sum_j C[i][j] atanh(1/x_j),

   and for the first n nonreal Gaussian primes a_i + b_i I,

       den arg(a_i + b_i I) = sum_j C[i][j] atan(1/x_j),

   each of which is exactly zero.  How many rows ca can PROVE zero is
   a measure of its integer-relation machinery: with the raised
   precision limits set below it proves all of them for the smaller
   sets and only some beyond about six terms, where the relation
   involves that many logarithms of large arguments.  A row reported
   as WRONG would mean a corrupt table. */
static int
check_machin_table(int gaussian, slong num, ca_ctx_t ctx)
{
    const mp_real_machin_struct * tab = _mp_real_machin_table(gaussian, num);
    ca_ptr series;
    ca_t x, y;
    fmpz * crow;
    fmpz_t xj;
    slong i, j, n, proved = 0;
    int ok = 1;

    n = tab->num;
    if (n != num)
        return 1;   /* no set of this size is tabulated */

    series = _ca_vec_init(n, ctx);
    crow = _fmpz_vec_init(n);
    ca_init(x, ctx);
    ca_init(y, ctx);
    fmpz_init(xj);

    /* the series of the set, atan(1/x_j) or atanh(1/x_j) */
    for (j = 0; j < n; j++)
    {
        _mp_real_machin_get_x(xj, tab, j);
        ca_one(series + j, ctx);
        ca_div_fmpz(series + j, series + j, xj, ctx);
        if (gaussian)
            ca_atan(series + j, series + j, ctx);
        else
            simple_ca_atanh(series + j, series + j, ctx);
    }

    for (i = 0; i < n; i++)
    {
        _mp_real_machin_get_c_row(crow, tab, i);

        ca_zero(x, ctx);
        for (j = 0; j < n; j++)
        {
            ca_mul_fmpz(y, series + j, crow + j, ctx);
            ca_add(x, x, y, ctx);
        }

        /* subtract den times the value the row represents */
        if (gaussian)
        {
            slong a = _mp_real_gaussian_primes[2 * i];
            slong b = _mp_real_gaussian_primes[2 * i + 1];
            ca_set_si(y, b, ctx);
            ca_div_si(y, y, a, ctx);
            ca_atan(y, y, ctx);
        }
        else
        {
            ca_set_ui(y, n_nth_prime(i + 1), ctx);
            ca_log(y, y, ctx);
        }
        ca_mul_ui(y, y, tab->den, ctx);
        ca_sub(x, x, y, ctx);

        {
            truth_t t = ca_check_is_zero(x, ctx);

            if (t == T_FALSE)
            {
                flint_printf("    row %wd is WRONG:   ", i);
                ca_print(x, ctx);
                flint_printf("\n");
                ok = -1;
            }
            else if (t == T_TRUE)
                proved++;
            else if (ok == 1)
                ok = 0;   /* ca could not decide; see below */
        }
    }

    /* ok = 1: every row proved zero; 0: some row left undecided (the
       formulas are verified numerically elsewhere, so this measures
       what ca can prove); -1: a row is provably nonzero */
    flint_printf("  %s, %2wd terms: %wd/%wd rows proved zero%s\n",
        gaussian ? "atan " : "atanh", n, proved, n,
        (ok == -1) ? "  (A FORMULA IS WRONG)" : "");

    _ca_vec_clear(series, n, ctx);
    _fmpz_vec_clear(crow, n);
    ca_clear(x, ctx);
    ca_clear(y, ctx);
    fmpz_clear(xj);
    return ok;
}

/* The same checks for a generic gr context (used with the lazy tower
   field, gr_ctx_init_tower_lazy, where logarithms of Gaussian rationals
   are decomposed over the logarithms of Gaussian primes, so that every
   row reduces to zero in the representation itself). */

/* atan(x) = (log(1 + i x) - log(1 - i x)) / (2 i), atanh(x) = log((1+x)/(1-x))/2 */
static int
gr_atan_via_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr i, t, u;
    int status = GR_SUCCESS;
    GR_TMP_INIT3(i, t, u, ctx);
    status |= gr_i(i, ctx);
    status |= gr_mul(t, i, x, ctx);
    status |= gr_add_ui(u, t, 1, ctx);
    status |= gr_log(u, u, ctx);
    status |= gr_neg(t, t, ctx);
    status |= gr_add_ui(t, t, 1, ctx);
    status |= gr_log(t, t, ctx);
    status |= gr_sub(u, u, t, ctx);
    status |= gr_mul_ui(t, i, 2, ctx);
    status |= gr_div(res, u, t, ctx);
    GR_TMP_CLEAR3(i, t, u, ctx);
    return status;
}

static int
gr_atanh_via_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t, u;
    int status = GR_SUCCESS;
    GR_TMP_INIT2(t, u, ctx);
    status |= gr_add_ui(t, x, 1, ctx);
    status |= gr_neg(u, x, ctx);
    status |= gr_add_ui(u, u, 1, ctx);
    status |= gr_div(res, t, u, ctx);
    status |= gr_log(res, res, ctx);
    status |= gr_div_ui(res, res, 2, ctx);
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

static int
check_machin_table_gr(int gaussian, slong num, gr_ctx_t ctx)
{
    const mp_real_machin_struct * tab = _mp_real_machin_table(gaussian, num);
    gr_ptr series, x, y;
    fmpz * crow;
    fmpz_t xj;
    slong i, j, n, proved = 0;
    int ok = 1;

    n = tab->num;
    if (n != num)
        return 1;

    series = gr_heap_init_vec(n, ctx);
    crow = _fmpz_vec_init(n);
    GR_TMP_INIT2(x, y, ctx);
    fmpz_init(xj);

    for (j = 0; j < n; j++)
    {
        gr_ptr sj = GR_ENTRY(series, j, ctx->sizeof_elem);
        _mp_real_machin_get_x(xj, tab, j);
        GR_MUST_SUCCEED(gr_one(sj, ctx));
        GR_MUST_SUCCEED(gr_div_fmpz(sj, sj, xj, ctx));
        if (gaussian)
            GR_MUST_SUCCEED(gr_atan_via_log(sj, sj, ctx));
        else
            GR_MUST_SUCCEED(gr_atanh_via_log(sj, sj, ctx));
    }

    for (i = 0; i < n; i++)
    {
        truth_t t;

        _mp_real_machin_get_c_row(crow, tab, i);
        GR_MUST_SUCCEED(gr_zero(x, ctx));
        for (j = 0; j < n; j++)
        {
            GR_MUST_SUCCEED(gr_mul_fmpz(y, GR_ENTRY(series, j, ctx->sizeof_elem), crow + j, ctx));
            GR_MUST_SUCCEED(gr_add(x, x, y, ctx));
        }

        if (gaussian)
        {
            slong a = _mp_real_gaussian_primes[2 * i];
            slong b = _mp_real_gaussian_primes[2 * i + 1];
            GR_MUST_SUCCEED(gr_set_si(y, b, ctx));
            GR_MUST_SUCCEED(gr_div_si(y, y, a, ctx));
            GR_MUST_SUCCEED(gr_atan_via_log(y, y, ctx));
        }
        else
        {
            GR_MUST_SUCCEED(gr_set_ui(y, n_nth_prime(i + 1), ctx));
            GR_MUST_SUCCEED(gr_log(y, y, ctx));
        }
        GR_MUST_SUCCEED(gr_mul_ui(y, y, tab->den, ctx));
        GR_MUST_SUCCEED(gr_sub(x, x, y, ctx));

        t = gr_is_zero(x, ctx);
        if (t == T_FALSE)
        {
            flint_printf("    row %wd is WRONG\n", i);
            ok = -1;
        }
        else if (t == T_TRUE)
            proved++;
        else if (ok == 1)
            ok = 0;
    }

    flint_printf("  %s, %2wd terms: %wd/%wd rows proved zero%s\n",
        gaussian ? "atan " : "atanh", n, proved, n,
        (ok == -1) ? "  (A FORMULA IS WRONG)" : "");

    gr_heap_clear_vec(series, n, ctx);
    _fmpz_vec_clear(crow, n);
    GR_TMP_CLEAR2(x, y, ctx);
    fmpz_clear(xj);
    return ok;
}

static void
run_gr(gr_ctx_t ctx, slong maxn)
{
    gr_ptr x, y, pi4;
    slong i, j, c, q;
    int g;
    slong num;

    GR_TMP_INIT3(x, y, pi4, ctx);
    GR_MUST_SUCCEED(gr_pi(pi4, ctx));
    GR_MUST_SUCCEED(gr_div_ui(pi4, pi4, 4, ctx));

    for (i = 0; i < NUM_FORMULAS; i++)
    {
        flint_printf("[");
        GR_MUST_SUCCEED(gr_zero(x, ctx));
        for (j = 0; j < 4; j++)
        {
            c = machin_formulas[i][j][0];
            q = machin_formulas[i][j][1];
            if (c != 0)
            {
                if (j != 0)
                    flint_printf(" + ");
                flint_printf("(%wd)*atan(1/%wd)", c, q);
                GR_MUST_SUCCEED(gr_one(y, ctx));
                GR_MUST_SUCCEED(gr_div_ui(y, y, q, ctx));
                GR_MUST_SUCCEED(gr_atan_via_log(y, y, ctx));
                GR_MUST_SUCCEED(gr_mul_si(y, y, c, ctx));
                GR_MUST_SUCCEED(gr_add(x, x, y, ctx));
            }
        }
        flint_printf(" - pi/4]   =   ");
        GR_MUST_SUCCEED(gr_sub(x, x, pi4, ctx));
        gr_println(x, ctx);
    }
    flint_printf("\n");

    for (i = 0; i < NUM_FORMULAS2; i++)
    {
        flint_printf("[");
        GR_MUST_SUCCEED(gr_zero(x, ctx));
        for (j = 0; j < 4; j++)
        {
            c = hyperbolic_machin_formulas[i][j][0];
            q = hyperbolic_machin_formulas[i][j][1];
            if (c != 0)
            {
                if (j != 0)
                    flint_printf(" + ");
                flint_printf("(%wd)*atanh(1/%wd)", c, q);
                GR_MUST_SUCCEED(gr_one(y, ctx));
                GR_MUST_SUCCEED(gr_div_ui(y, y, q, ctx));
                GR_MUST_SUCCEED(gr_atanh_via_log(y, y, ctx));
                GR_MUST_SUCCEED(gr_mul_si(y, y, c, ctx));
                GR_MUST_SUCCEED(gr_add(x, x, y, ctx));
            }
        }
        flint_printf(" - log(%wd)]   =   ", hyperbolic_logs[i]);
        GR_MUST_SUCCEED(gr_set_ui(y, hyperbolic_logs[i], ctx));
        GR_MUST_SUCCEED(gr_log(y, y, ctx));
        GR_MUST_SUCCEED(gr_sub(x, x, y, ctx));
        gr_println(x, ctx);
    }

    flint_printf("\nFLINT's tabulated Machin-type sets"
        " (_mp_real_machin_table), sum - den * value:\n");
    for (g = 0; g < 2; g++)
        for (num = 2; num <= maxn; num++)
            check_machin_table_gr(g, num, ctx);

    GR_TMP_CLEAR3(x, y, pi4, ctx);
}

int main(int argc, char *argv[])

{
    ca_ctx_t ctx;
    ca_t x, y, pi4;
    slong i, j, c, q;

    if (argc > 1 && strcmp(argv[1], "-tower") == 0)
    {
        gr_ctx_t QQ, K;
        slong maxn = (argc > 2) ? atol(argv[2]) : 48;
        TIMEIT_ONCE_START;
        gr_ctx_init_fmpq(QQ);
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        run_gr(K, maxn);
        gr_tower_lazy_ctx_stats(K);
        gr_ctx_clear(K);
        gr_ctx_clear(QQ);
        flint_printf("\n");
        TIMEIT_ONCE_STOP;
        print_memory_usage();
        flint_cleanup();
        return 0;
    }

    TIMEIT_ONCE_START;

    ca_ctx_init(ctx);
    ca_init(x, ctx);
    ca_init(y, ctx);
    ca_init(pi4, ctx);

    ca_pi(pi4, ctx);
    ca_div_ui(pi4, pi4, 4, ctx);

    for (i = 0; i < NUM_FORMULAS; i++)
    {
        flint_printf("[");
        ca_zero(x, ctx);
        for (j = 0; j < 4; j++)
        {
            c = machin_formulas[i][j][0];
            q = machin_formulas[i][j][1];

            if (c != 0)
            {
                if (j != 0)
                    flint_printf(" + ");
                flint_printf("(%wd)*atan(1/%wd)", c, q);
                simple_ca_atan_p_q(y, 1, q, ctx);
                ca_mul_si(y, y, c, ctx);
                ca_add(x, x, y, ctx);
            }
        }

        flint_printf(" - pi/4]   =   ");
        ca_sub(x, x, pi4, ctx);

        ca_print(x, ctx);
        flint_printf("\n");
    }

    flint_printf("\n");

    for (i = 0; i < NUM_FORMULAS2; i++)
    {
        flint_printf("[");
        ca_zero(x, ctx);
        for (j = 0; j < 4; j++)
        {
            c = hyperbolic_machin_formulas[i][j][0];
            q = hyperbolic_machin_formulas[i][j][1];

            if (c != 0)
            {
                if (j != 0)
                    flint_printf(" + ");
                flint_printf("(%wd)*atanh(1/%wd)", c, q);
                simple_ca_atanh_p_q(y, 1, q, ctx);
                ca_mul_si(y, y, c, ctx);
                ca_add(x, x, y, ctx);
            }
        }

        flint_printf(" - log(%wd)]   =   ", hyperbolic_logs[i]);
        ca_set_ui(y, hyperbolic_logs[i], ctx);
        ca_log(y, y, ctx);
        ca_sub(x, x, y, ctx);

        ca_print(x, ctx);
        flint_printf("\n");
    }

    /* the tabulated sets: every size up to maxn */
    {
        slong maxn = (argc > 1) ? atol(argv[1]) : 8;
        slong num;
        int g;
        ca_ctx_t ctx2;

        /* deciding these needs integer relation detection among many
           logarithms, for which the default LLL and evaluation
           precisions are not enough beyond about five terms; the
           options must be set on a fresh context, before any
           extension objects have been cached */
        ca_ctx_init(ctx2);
        ctx2->options[CA_OPT_LLL_PREC] = 4096;
        ctx2->options[CA_OPT_PREC_LIMIT] = 65536;

        flint_printf("\nFLINT's tabulated Machin-type sets"
            " (_mp_real_machin_table), sum - den * value:\n");
        for (g = 0; g < 2; g++)
            for (num = 2; num <= maxn; num++)
                check_machin_table(g, num, ctx2);

        ca_ctx_clear(ctx2);
    }

    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_clear(pi4, ctx); 
    ca_ctx_clear(ctx);

    flint_printf("\n");
    TIMEIT_ONCE_STOP;
    print_memory_usage();

    flint_cleanup();
    return 0;
}
