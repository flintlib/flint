/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/acb_modular.h>
#include <flint/fmpz_poly.h>
#include <flint/profiler.h>

int main(int argc, char *argv[])
{
    fmpz_poly_t res;
    slong i, num_threads;
    slong D;

    if (argc < 2)
    {
        flint_printf("usage: build/examples/class_poly D [-threads n]\n");
        return 1;
    }

    D = atol(argv[1]);

    num_threads = 1;

    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
    }

    flint_set_num_threads(num_threads);

    fmpz_poly_init(res);

    TIMEIT_ONCE_START
    acb_modular_hilbert_class_poly(res, D);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    if (FLINT_ABS(D) <= 100)
    {
        fmpz_poly_print(res);
        flint_printf("\n");
    }
    else
    {
        flint_printf("degree = %wd, bits = %wd\n",
            fmpz_poly_degree(res), fmpz_poly_max_bits(res));
    }

    fmpz_poly_clear(res);
    flint_cleanup_master();
    return 0;
}

