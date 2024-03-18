/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <stdlib.h>
#include "arb.h"
#include "profiler.h"

int main(int argc, char *argv[])
{
    arb_t x;
    slong i, prec, digits, condense, num_threads;
    int constant;

    if (argc < 2)
    {
        flint_printf("usage: build/examples/pi digits [-condense n] [-threads n] [-constant name]\n");
        flint_printf("compute digits of pi (or other constants)\n");
        flint_printf("supported: pi, e, log2, euler, catalan, khinchin, zeta3, reciprocal_fibonacci\n");
        return 1;
    }

    digits = atol(argv[1]);

    num_threads = 1;
    condense = 20;
    constant = 0;

    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-condense"))
            condense = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-constant"))
        {
            if (!strcmp(argv[i+1], "pi"))
                constant = 0;
            else if (!strcmp(argv[i+1], "e"))
                constant = 1;
            else if (!strcmp(argv[i+1], "log2"))
                constant = 2;
            else if (!strcmp(argv[i+1], "euler"))
                constant = 3;
            else if (!strcmp(argv[i+1], "catalan"))
                constant = 4;
            else if (!strcmp(argv[i+1], "khinchin"))
                constant = 5;
            else if (!strcmp(argv[i+1], "zeta3"))
                constant = 6;
            else if (!strcmp(argv[i+1], "reciprocal_fibonacci"))
                constant = 7;
            else
            {
                flint_printf("unknown constant\n");
                flint_abort();
            }
        }
    }

    flint_set_num_threads(num_threads);

    arb_init(x);

    prec = digits * 3.3219280948873623 + 5;

    flint_printf("precision = %wd bits... ", prec);
    fflush(stdout);

    TIMEIT_ONCE_START
    if (constant == 0)
        arb_const_pi(x, prec);
    else if (constant == 1)
        arb_const_e(x, prec);
    else if (constant == 2)
        arb_const_log2(x, prec);
    else if (constant == 3)
        arb_const_euler(x, prec);
    else if (constant == 4)
        arb_const_catalan(x, prec);
    else if (constant == 5)
        arb_const_khinchin(x, prec);
    else if (constant == 6)
        arb_zeta_ui(x, 3, prec);
    else if (constant == 7)
        arb_const_reciprocal_fibonacci(x, prec);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    arb_printn(x, digits, ARB_STR_CONDENSE * condense);

    flint_printf("\n");

    arb_clear(x);
    flint_cleanup_master();
    return 0;
}

