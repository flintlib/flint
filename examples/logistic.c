/* This file is public domain. Author: Fredrik Johansson. */

#include "arb.h"
#include "flint/profiler.h"

int main(int argc, char *argv[])
{
    slong prec, goal, digits;
    const char * xs;
    const char * rs;
    arb_t x, r, s, t;
    slong i, n;

    if (argc < 2)
    {
        flint_printf("nth iterate of the logistic map x_{n+1} = r x_n (1-x_n)\n");
        flint_printf("usage: build/examples/logistic n [x_0=0.5] [r=3.75] [digits=10]\n");
        return 1;
    }

    n = atol(argv[1]);
    n = FLINT_MAX(n, 0);
    xs = (argc >= 3) ? argv[2] : "0.5";
    rs = (argc >= 4) ? argv[3] : "3.75";
    digits = (argc >= 5) ? atol(argv[4]) : 10;
    if (digits < 1)
        digits = 1;

    goal = digits * 3.3219280948873623 + 3;

    arb_init(x);
    arb_init(r);
    arb_init(s);
    arb_init(t);

    TIMEIT_ONCE_START
    for (prec = 64; ; prec *= 2)
    {
        flint_printf("Trying prec=%wd bits...", prec);
        fflush(stdout);

        arb_set_str(x, xs, prec);
        arb_set_str(r, rs, prec);

        if (!arb_is_finite(x) || !arb_is_finite(r))
        {
            flint_printf("unable to parse input string!\n");
            flint_abort();
        }

        for (i = 0; i < n; i++)
        {
            arb_sub_ui(t, x, 1, prec);
            arb_neg(t, t);
            arb_mul(x, x, t, prec);
            arb_mul(x, x, r, prec);

            if (arb_rel_accuracy_bits(x) < goal)
            {
                flint_printf("ran out of accuracy at step %wd\n", i);
                break;
            }
        }

        if (i == n)
        {
            flint_printf("success!\n");
            break;
        }
    }
    TIMEIT_ONCE_STOP

    flint_printf("x_%wd = ", n);
    arb_printn(x, digits, 0);
    flint_printf("\n");

    arb_clear(x);
    arb_clear(r);
    arb_clear(s);
    arb_clear(t);

    flint_cleanup();
    return 0;
}

