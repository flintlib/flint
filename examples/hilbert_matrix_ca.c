/* This file is public domain. Author: Fredrik Johansson. */

#include <stdlib.h>
#include <string.h>
#include <flint/profiler.h>
#include <flint/fmpq_mat.h>
#include <flint/calcium.h>
#include <flint/ca.h>
#include <flint/ca_vec.h>
#include <flint/ca_mat.h>
#include <flint/qqbar.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/gr_poly.h>
#include <flint/gr_vec.h>
#include <flint/fmpz_vec.h>
#include <flint/gr_tower.h>
#include <flint/gr_tower_lazy.h>

int main(int argc, char *argv[])
{
    slong n, i;
    int qqbar, vieta, novieta, tower;

    if (argc < 2)
    {
        flint_printf("usage: hilbert_matrix [-qqbar] [-tower] [-vieta | -novieta] n\n");
        return 1;
    }

    qqbar = 0;
    tower = 0;
    vieta = 0;
    novieta = 0;
    n = 0;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-qqbar"))
        {
            qqbar = 1;
        }
        else if (!strcmp(argv[i], "-tower"))
        {
            tower = 1;
        }
        else if (!strcmp(argv[i], "-vieta"))
        {
            vieta = 1;
        }
        else if (!strcmp(argv[i], "-novieta"))
        {
            novieta = 1;
        }
        else
        {
            n = atol(argv[i]);
            if (n < 0 || n > 100)
                flint_abort();
        }
    }

    TIMEIT_ONCE_START;

    if (qqbar)
    {
        fmpq_mat_t mat;
        qqbar_ptr eig;
        qqbar_t trace, det;

        fmpq_mat_init(mat, n, n);
        qqbar_init(trace);
        qqbar_init(det);
        eig = _qqbar_vec_init(n);

        fmpq_mat_hilbert_matrix(mat);
        qqbar_eigenvalues_fmpq_mat(eig, mat, 0);

        flint_printf("Trace:\n");
        qqbar_zero(trace);
        for (i = 0; i < n; i++)
        {
            qqbar_add(trace, trace, eig + i);
            flint_printf("%wd/%wd: degree %wd\n", i, n, qqbar_degree(trace));
        }
        qqbar_print(trace);
        flint_printf("\n");

        flint_printf("Determinant:\n");
        qqbar_one(det);
        for (i = 0; i < n; i++)
        {
            qqbar_mul(det, det, eig + i);
            flint_printf("%wd/%wd: degree %wd\n", i, n, qqbar_degree(det));
        }
        qqbar_print(det);
        flint_printf("\n");

        fmpq_mat_clear(mat);
        qqbar_clear(trace);
        qqbar_clear(det);
        _qqbar_vec_clear(eig, n);
    }
    else if (tower)
    {
        /* the lazy tower field: the eigenvalues of the irreducible
           factors of the characteristic polynomial generate splitting
           towers, in which symmetric functions of the roots reduce to
           the coefficients */
        gr_ctx_t QQ, K;
        gr_mat_t mat;
        gr_poly_t cp;
        gr_vec_t eig;
        fmpz_vec_t mul;
        gr_ptr trace, det, t;
        fmpq_mat_t hmat;

        gr_ctx_init_fmpq(QQ);
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);

        fmpq_mat_init(hmat, n, n);
        fmpq_mat_hilbert_matrix(hmat);
        gr_mat_init(mat, n, n, QQ);
        GR_MUST_SUCCEED(gr_mat_set_fmpq_mat(mat, hmat, QQ));

        gr_poly_init(cp, QQ);
        gr_vec_init(eig, 0, K);
        fmpz_vec_init(mul, 0);
        trace = gr_heap_init(K);
        det = gr_heap_init(K);
        t = gr_heap_init(K);

        GR_MUST_SUCCEED(gr_mat_charpoly(cp, mat, QQ));

        {
            gr_poly_t cpK;
            gr_poly_init(cpK, K);
            GR_MUST_SUCCEED(gr_poly_set_gr_poly_other(cpK, cp, QQ, K));
            GR_MUST_SUCCEED(gr_poly_roots(eig, mul, cpK, 0, K));
            gr_poly_clear(cpK, K);
        }

        {
            gr_ptr tr, dt;
            tr = gr_heap_init(QQ);
            dt = gr_heap_init(QQ);
            GR_MUST_SUCCEED(gr_mat_trace(tr, mat, QQ));
            GR_MUST_SUCCEED(gr_mat_det(dt, mat, QQ));
            GR_MUST_SUCCEED(gr_set_other(trace, tr, QQ, K));
            GR_MUST_SUCCEED(gr_set_other(det, dt, QQ, K));
            gr_heap_clear(tr, QQ);
            gr_heap_clear(dt, QQ);
        }

        flint_printf("Trace:\n");
        GR_MUST_SUCCEED(gr_zero(t, K));
        for (i = 0; i < eig->length; i++)
        {
            slong j;
            for (j = 0; j < fmpz_get_si(mul->entries + i); j++)
                GR_MUST_SUCCEED(gr_add(t, t, gr_vec_entry_ptr(eig, i, K), K));
        }
        gr_println(trace, K);
        gr_println(t, K);
        flint_printf("Equal: "); truth_print(gr_equal(trace, t, K)); flint_printf("\n\n");

        flint_printf("Det:\n");
        GR_MUST_SUCCEED(gr_one(t, K));
        for (i = 0; i < eig->length; i++)
        {
            slong j;
            for (j = 0; j < fmpz_get_si(mul->entries + i); j++)
                GR_MUST_SUCCEED(gr_mul(t, t, gr_vec_entry_ptr(eig, i, K), K));
        }
        gr_println(det, K);
        gr_println(t, K);
        flint_printf("Equal: "); truth_print(gr_equal(det, t, K)); flint_printf("\n\n");

        gr_tower_lazy_ctx_stats(K);

        gr_heap_clear(trace, K);
        gr_heap_clear(det, K);
        gr_heap_clear(t, K);
        gr_vec_clear(eig, K);
        fmpz_vec_clear(mul);
        gr_poly_clear(cp, QQ);
        gr_mat_clear(mat, QQ);
        fmpq_mat_clear(hmat);
        gr_ctx_clear(K);
        gr_ctx_clear(QQ);
    }
    else
    {
        ca_ctx_t ctx;
        ca_mat_t mat;
        ca_vec_t eig;
        ca_t trace, det, t;
        ulong * mul;

        ca_ctx_init(ctx);

        /* Verification requires high-degree algebraics. */
        ctx->options[CA_OPT_QQBAR_DEG_LIMIT] = 10000;

        if (vieta)
            ctx->options[CA_OPT_VIETA_LIMIT] = n;
        if (novieta)
            ctx->options[CA_OPT_VIETA_LIMIT] = 0;

        ca_mat_init(mat, n, n, ctx);
        ca_vec_init(eig, 0, ctx);
        mul = flint_malloc(sizeof(ulong) * n);
        ca_init(trace, ctx);
        ca_init(det, ctx);
        ca_init(t, ctx);

        ca_mat_hilbert(mat, ctx);

        ca_mat_eigenvalues(eig, mul, mat, ctx);

        ca_mat_trace(trace, mat, ctx);
        ca_mat_det(det, mat, ctx);

        /* note: in general, we should use the multiplicities, but
           we happen to know that the eigenvalues are simple here */

        flint_printf("Trace:\n");
        ca_zero(t, ctx);
        for (i = 0; i < n; i++)
            ca_add(t, t, ca_vec_entry(eig, i), ctx);

        ca_print(trace, ctx); flint_printf("\n");
        ca_print(t, ctx); flint_printf("\n");
        flint_printf("Equal: "); truth_print(ca_check_equal(trace, t, ctx)); flint_printf("\n\n");

        flint_printf("Det:\n");
        ca_one(t, ctx);
        for (i = 0; i < n; i++)
            ca_mul(t, t, ca_vec_entry(eig, i), ctx);

        ca_print(det, ctx); flint_printf("\n");
        ca_print(t, ctx); flint_printf("\n");
        flint_printf("Equal: "); truth_print(ca_check_equal(det, t, ctx)); flint_printf("\n\n");

        ca_mat_clear(mat, ctx);
        ca_vec_clear(eig, ctx);
        flint_free(mul);
        ca_clear(trace, ctx);
        ca_clear(det, ctx);
        ca_clear(t, ctx);
        ca_ctx_clear(ctx);
    }

    flint_printf("\n");
    TIMEIT_ONCE_STOP;
    print_memory_usage();

    flint_cleanup();
    return 0;
}
