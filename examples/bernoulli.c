/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <stdlib.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/bernoulli.h>
#include <flint/profiler.h>

int main(int argc, char *argv[])
{
    fmpq_t x;
    slong i, num_threads;
    ulong n;

    if (argc < 2)
    {
        flint_printf("usage: build/examples/bernoulli n [-threads n]\n");
        return 1;
    }

    n = atol(argv[1]);

    num_threads = 1;

    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
    }

    flint_set_num_threads(num_threads);

    fmpq_init(x);

    TIMEIT_ONCE_START
    bernoulli_fmpq_ui(x, n);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    if (n <= 100)
    {
        fmpq_print(x);
        flint_printf("\n");
    }

    fmpq_clear(x);
    flint_cleanup_master();
    return 0;
}

