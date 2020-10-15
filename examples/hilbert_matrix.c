/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "flint/fmpq_mat.h"
#include "ca.h"
#include "ca_vec.h"
#include "ca_mat.h"
#include "qqbar.h"

int main(int argc, char *argv[])
{
    slong n, i;
    int qqbar, vieta, novieta;

    if (argc < 2)
    {
        flint_printf("usage: hilbert_matrix [-qqbar] [-vieta | -novieta] n\n");
        return 1;
    }

    qqbar = 0;
    vieta = 0;
    novieta = 0;
    n = 0;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-qqbar"))
        {
            qqbar = 1;
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

    TIMEIT_ONCE_START

    if (qqbar)
    {
        fmpq_mat_t mat;
        qqbar_ptr eig;
        qqbar_t trace, det;

        fmpq_mat_init(mat, n, n);
        qqbar_init(trace);
        qqbar_init(det);
        eig = qqbar_vec_init(n);

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
        qqbar_vec_clear(eig, n);
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
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}
