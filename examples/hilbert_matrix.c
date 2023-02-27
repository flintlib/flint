/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "arb_mat.h"
#include "acb_mat.h"
#include "flint/profiler.h"

int main(int argc, char *argv[])
{
    arb_mat_t A;
    arb_t det;
    slong i, prec, n;
    int eig;

    if (argc < 2 || (argc == 3 && strcmp(argv[1], "-eig")))
    {
        flint_printf("usage: build/examples/hilbert_matrix [-eig] n\n");
        return 1;
    }

    if (argc == 2)
        eig = 0;
    else
        eig = 1;

    n = atol(argv[argc-1]);

    if (eig && (n == 0))
    {
        flint_printf("smallest eigenvalue: undefined for n == 0\n");
        return 1;
    }

    arb_mat_init(A, n, n);
    arb_init(det);

    TIMEIT_ONCE_START

    for (prec = 20; ; prec *= 2)
    {
        arb_mat_hilbert(A, prec);
        flint_printf("prec=%wd: ", prec);

        if (eig == 0)
        {
            arb_mat_det(det, A, prec);
        }
        else
        {
            acb_mat_t C, R;
            acb_ptr E;

            acb_mat_init(R, n, n);
            acb_mat_init(C, n, n);
            E = _acb_vec_init(n);

            acb_mat_set_arb_mat(C, A);
            acb_mat_approx_eig_qr(E, NULL, R, C, NULL, 0, prec);
            if (acb_mat_eig_simple(E, NULL, NULL, C, E, R, prec))
            {
                /* A is symmetric so the eigenvalues are real */
                for (i = 0; i < n; i++)
                    arb_zero(acb_imagref(E + i));
                /* With isolated eigenvalues, sorting midpoints gives the
                   right result. */
                _acb_vec_sort_pretty(E, n);
                acb_get_real(det, E + 0);
            }
            else
            {
                arb_indeterminate(det);
            }

            acb_mat_clear(R);
            acb_mat_clear(C);
            _acb_vec_clear(E, n);
        }

        arb_printn(det, 10, 0);
        flint_printf("\n");

        if (!arb_contains_zero(det))
        {
            flint_printf("success!\n");
            break;
        }
    }

    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    arb_mat_clear(A);
    arb_clear(det);
    flint_cleanup();
    return 0;
}
