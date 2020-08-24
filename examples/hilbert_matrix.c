/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "flint/fmpq_mat.h"
#include "qqbar.h"

int main(int argc, char *argv[])
{
    fmpq_mat_t mat;
    qqbar_ptr eig;
    qqbar_t trace, det;
    slong n, i;

    if (argc < 2)
    {
        flint_printf("usage: hilbert_matrix n\n");
        return 1;
    }

    n = atol(argv[1]);
    if (n < 0 || n > 100)
        flint_abort();

    TIMEIT_ONCE_START

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

    flint_printf("\n");
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}
