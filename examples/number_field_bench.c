/* This file is public domain. Author: Fredrik Johansson. */

/*
    Arithmetic and linear algebra over number fields Q(zeta_N) and
    Q(sqrt(D)) in three gr contexts: Calcium (gr_ctx_init_complex_ca), the
    lazy tower field (gr_ctx_init_tower_lazy) and nf (gr_ctx_init_nf,
    which knows the field in advance). Elements are dense random
    combinations of the powers of the generator with small coefficients.

    usage: number_field_bench [-ops] [-linalg n] N ...

    N > 0 selects Q(zeta_N), N < 0 selects Q(sqrt(N)). With -ops (the
    default): products of elements, polynomial products (length 20),
    matrix products (10 x 10), determinants and linear systems (6 x 6).
    With -linalg n: determinants and linear systems of size n with each
    algorithm (LU, fraction-free LU, Berkowitz; the default choice of the
    context last).
*/

#include <stdlib.h>
#include <string.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_poly.h>
#include <flint/gr_mat.h>
#include <flint/gr_special.h>
#include <flint/gr_tower.h>
#include <flint/gr_tower_lazy.h>

static const char * ctx_names[] = {"ca", "tower", "nf"};

static void
init_ctx(gr_ctx_t ctx, gr_ctx_t QQ, int which, slong N)
{
    if (which == 0)
        gr_ctx_init_complex_ca(ctx);
    else if (which == 1)
        gr_ctx_init_tower_lazy(ctx, QQ, GR_TOWER_MERGE_EXPRESS);
    else
    {
        fmpz_poly_t m;
        fmpq_poly_t mq;
        fmpz_poly_init(m);
        fmpq_poly_init(mq);
        if (N > 0)
            fmpz_poly_cyclotomic(m, N);
        else
        {
            fmpz_poly_set_coeff_si(m, 2, 1);
            fmpz_poly_set_coeff_si(m, 0, -N);
        }
        fmpq_poly_set_fmpz_poly(mq, m);
        gr_ctx_init_nf(ctx, mq);
        fmpz_poly_clear(m);
        fmpq_poly_clear(mq);
    }
}

/* the generator: exp(2 pi i / N) or sqrt(N) */
static void
generator(gr_ptr z, slong N, int which, gr_ctx_t ctx)
{
    if (which == 2)
    {
        GR_MUST_SUCCEED(gr_gen(z, ctx));
    }
    else if (N > 0)
    {
        gr_ptr t = gr_heap_init(ctx);
        GR_MUST_SUCCEED(gr_pi(z, ctx));
        GR_MUST_SUCCEED(gr_i(t, ctx));
        GR_MUST_SUCCEED(gr_mul(z, z, t, ctx));
        GR_MUST_SUCCEED(gr_mul_ui(z, z, 2, ctx));
        GR_MUST_SUCCEED(gr_div_si(z, z, N, ctx));
        GR_MUST_SUCCEED(gr_exp(z, z, ctx));
        gr_heap_clear(t, ctx);
    }
    else
    {
        GR_MUST_SUCCEED(gr_set_si(z, N, ctx));
        GR_MUST_SUCCEED(gr_sqrt(z, z, ctx));
    }
}

static void
rand_elem(gr_ptr x, gr_srcptr z, slong d, flint_rand_t state, gr_ctx_t ctx)
{
    slong j;
    GR_MUST_SUCCEED(gr_zero(x, ctx));
    for (j = d - 1; j >= 0; j--)
    {
        GR_MUST_SUCCEED(gr_mul(x, x, z, ctx));
        GR_MUST_SUCCEED(gr_add_si(x, x, (slong) n_randint(state, 7) - 3, ctx));
    }
}

static void
rand_mat(gr_mat_t A, gr_srcptr z, slong d, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    for (i = 0; i < gr_mat_nrows(A, ctx); i++)
        for (j = 0; j < gr_mat_ncols(A, ctx); j++)
            rand_elem(gr_mat_entry_ptr(A, i, j, ctx), z, d, state, ctx);
}

#define NUM_OPS 5
static const char * op_names[NUM_OPS] = {"mul", "poly_mul", "mat_mul", "det", "solve"};

#define NUM_ALGS 6
static const char * alg_names[NUM_ALGS] = {"det_lu", "det_fflu", "det_berkowitz", "solve_lu", "solve_fflu", "det/solve"};

/* time of operation k with the context `which`, in seconds (-1 if it failed) */
static double
run(slong N, int linalg, slong n, slong k, int which)
{
    gr_ctx_t ctx, QQ;
    gr_ptr z, x, y, r;
    flint_rand_t state;
    slong d = (N > 0) ? (slong) n_euler_phi(N) : 2, i;
    timeit_t timer;
    slong reps;
    int status = GR_SUCCESS;
    double t;

    flint_rand_init(state);
    gr_ctx_init_fmpq(QQ);
    init_ctx(ctx, QQ, which, N);
    z = gr_heap_init(ctx);
    x = gr_heap_init(ctx);
    y = gr_heap_init(ctx);
    r = gr_heap_init(ctx);
    generator(z, N, which, ctx);

    if (!linalg && k == 0)
    {
        gr_ptr x0 = gr_heap_init(ctx);
        rand_elem(x0, z, d, state, ctx);
        TIMEIT_REPEAT(timer, reps)
        {
            flint_rand_t st2;
            flint_rand_init(st2);
            status |= gr_set(x, x0, ctx);
            for (i = 0; i < 100; i++)
            {
                rand_elem(y, z, d, st2, ctx);
                status |= gr_mul(r, x, y, ctx);
                status |= gr_add(x, x, r, ctx);
                status |= gr_div_ui(x, x, 3, ctx);
            }
            flint_rand_clear(st2);
        }
        TIMEIT_END_REPEAT(timer, reps);
        gr_heap_clear(x0, ctx);
    }
    else if (!linalg && k == 1)
    {
        gr_poly_t f, g, h;
        slong len = 20;
        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(h, ctx);
        gr_poly_fit_length(f, len, ctx);
        gr_poly_fit_length(g, len, ctx);
        for (i = 0; i < len; i++)
        {
            rand_elem(gr_poly_coeff_ptr(f, i, ctx), z, d, state, ctx);
            rand_elem(gr_poly_coeff_ptr(g, i, ctx), z, d, state, ctx);
        }
        _gr_poly_set_length(f, len, ctx);
        _gr_poly_set_length(g, len, ctx);
        TIMEIT_REPEAT(timer, reps)
            status |= gr_poly_mul(h, f, g, ctx);
        TIMEIT_END_REPEAT(timer, reps);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(h, ctx);
    }
    else
    {
        gr_mat_t A, B, X;
        slong m = linalg ? n : ((k == 2) ? 10 : 6);
        gr_mat_init(A, m, m, ctx);
        gr_mat_init(B, m, (!linalg && k == 2) ? m : 1, ctx);
        gr_mat_init(X, m, gr_mat_ncols(B, ctx), ctx);
        rand_mat(A, z, d, state, ctx);
        rand_mat(B, z, d, state, ctx);

        TIMEIT_REPEAT(timer, reps)
        if (!linalg)
        {
            if (k == 2)
                status |= gr_mat_mul(X, A, B, ctx);
            else if (k == 3)
                status |= gr_mat_det(r, A, ctx);
            else
                status |= gr_mat_nonsingular_solve(X, A, B, ctx);
        }
        else
        {
            if (k == 0)
                status |= gr_mat_det_lu(r, A, ctx);
            else if (k == 1)
                status |= gr_mat_det_fflu(r, A, ctx);
            else if (k == 2)
                status |= gr_mat_det_berkowitz(r, A, ctx);
            else if (k == 3)
                status |= gr_mat_nonsingular_solve_lu(X, A, B, ctx);
            else if (k == 4)
                status |= gr_mat_nonsingular_solve_fflu(X, A, B, ctx);
            else
            {
                status |= gr_mat_det(r, A, ctx);
                status |= gr_mat_nonsingular_solve(X, A, B, ctx);
            }
        }
        TIMEIT_END_REPEAT(timer, reps);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(X, ctx);
    }

    t = (status == GR_SUCCESS) ? timer->cpu * 0.001 / reps : -1.0;

    gr_heap_clear(z, ctx);
    gr_heap_clear(x, ctx);
    gr_heap_clear(y, ctx);
    gr_heap_clear(r, ctx);
    gr_ctx_clear(ctx);
    gr_ctx_clear(QQ);
    flint_rand_clear(state);
    return t;
}

int main(int argc, char *argv[])
{
    int linalg = 0, which, have_field = 0;
    slong n = 8, i, k;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-ops"))
            linalg = 0;
        else if (!strcmp(argv[i], "-linalg") && i + 1 < argc)
        {
            linalg = 1;
            n = atol(argv[++i]);
        }
        else
        {
            slong N = atol(argv[i]);

            if (N == 0 || N == 1 || N == 2)
            {
                flint_printf("usage: number_field_bench [-ops] [-linalg n] N ...  (N > 2: Q(zeta_N), N < 0: Q(sqrt(N)))\n");
                return 1;
            }

            have_field = 1;
            if (N > 0)
                flint_printf("Q(zeta_%wd), degree %wu", N, n_euler_phi(N));
            else
                flint_printf("Q(sqrt(%wd)), degree 2", N);
            if (linalg)
                flint_printf(", %wd x %wd matrices", n, n);
            flint_printf(" (times in ms)\n%-16s %10s %10s %10s\n", "", ctx_names[0], ctx_names[1], ctx_names[2]);

            for (k = 0; k < (linalg ? NUM_ALGS : NUM_OPS); k++)
            {
                flint_printf("%-16s", linalg ? alg_names[k] : op_names[k]);
                for (which = 0; which < 3; which++)
                {
                    /* (Berkowitz beyond 16 x 16 is too slow to be interesting) */
                    if (linalg && k == 2 && n > 16)
                        flint_printf(" %10s", "-");
                    else
                    {
                        double t = run(N, linalg, n, k, which);
                        if (t < 0)
                            flint_printf(" %10s", "failed");
                        else
                            flint_printf(" %10.4g", 1000 * t);
                    }
                    fflush(stdout);
                }
                flint_printf("\n");
            }
            flint_printf("\n");
        }
    }

    if (!have_field)
    {
        flint_printf("usage: number_field_bench [-ops] [-linalg n] N ...  (N > 2: Q(zeta_N), N < 0: Q(sqrt(N)))\n");
        return 1;
    }

    flint_cleanup_master();
    return 0;
}
