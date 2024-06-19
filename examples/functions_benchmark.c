/* This file is public domain. Author: Fredrik Johansson. */

#include <flint/profiler.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#include <flint/acb.h>
#include <flint/acb_modular.h>
#include <flint/partitions.h>
#include <flint/bernoulli.h>

#define TIMEIT_ONCE_STOP_VALUES(tcpu, twall) \
        } while (0); \
        timeit_stop(__timer); \
        (tcpu) = __timer->cpu*0.001; \
        (twall) = __timer->wall*0.001; \
    } while (0);

void doit(arb_t res, const arb_t x, slong n, int function, slong prec)
{
    if (function == 0)
        arb_const_pi(res, prec);
    else if (function == 1)
        arb_const_euler(res, prec);
    else if (function == 2)
        arb_exp(res, x, prec);
    else if (function == 3)
        arb_log(res, x, prec);
    else if (function == 4)
        arb_sin(res, x, prec);
    else if (function == 5)
        arb_atan(res, x, prec);
    else if (function == 6)
        arb_hypgeom_erf(res, x, prec);
    else if (function == 7)
        arb_gamma(res, x, prec);
    else if (function == 8)
        arb_zeta(res, x, prec);
    else if (function == 9)
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_imagref(t), x);
        acb_modular_eta(t, t, prec);
        acb_clear(t);
    }
    else if (function == 10)
    {
        fmpq_t t;
        fmpq_init(t);
        bernoulli_fmpq_ui(t, n);
        fmpq_clear(t);
    }
    else if (function == 11)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_ui(t, n);
        fmpz_mul(t, t, t);
        partitions_fmpz_fmpz(t, t, 0);
        fmpz_clear(t);
    }
}

char * description[] = {
    "const_pi, n digits",
    "const_euler, n digits",
    "exp(x), n digits",
    "log(x), n digits",
    "sin(x), n digits",
    "atan(x), n digits",
    "erf(x), n digits",
    "gamma(x), n digits",
    "zeta(x), n digits",
    "eta(ix), n digits",
    "bernoulli(n)",
    "partitions(n^2)",
};

slong limit[] = {
    1e8,
    1e7,
    1e7,
    1e7,
    1e7,
    1e7,
    1e6,
    1e5,
    1e4,
    1e6,
    1e6,
    1e7,
};

int main(void)
{
    arb_t x, y, res;
    slong n, prec;
    int function;
    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    arb_init(x);
    arb_init(y);
    arb_init(res);

    for (function = 0; function <= 11; function++)
    {
        printf("\n%s\n", description[function]);
        printf("%12s%24s%24s\n", "", "threads = 1   ", "threads = 8   ");
        printf("%12s%12s%12s%12s%12s\n", "n", "first", "repeat", "first", "repeat");

        for (n = 10; n <= limit[function]; n *= 10)
        {
            prec = n * 3.323 + 1;

            flint_printf("%12wd", n);
            fflush(stdout);

            arb_sqrt_ui(x, 2, prec);
            arb_sub_ui(x, x, 1, prec);

            flint_set_num_threads(1);

            flint_cleanup();
            TIMEIT_ONCE_START
            doit(res, x, n, function, prec);
            TIMEIT_ONCE_STOP_VALUES(tcpu, twall);
            printf("%12.3g", twall);
            fflush(stdout);

            TIMEIT_START
            doit(res, x, n, function, prec);
            TIMEIT_STOP_VALUES(tcpu, twall);
            printf("%12.3g", twall);
            fflush(stdout);

            flint_set_num_threads(8);

            flint_cleanup();
            TIMEIT_ONCE_START
            doit(res, x, n, function, prec);
            TIMEIT_ONCE_STOP_VALUES(tcpu, twall);
            printf("%12.3g", twall);
            fflush(stdout);

            TIMEIT_START
            doit(res, x, n, function, prec);
            TIMEIT_STOP_VALUES(tcpu, twall);
            printf("%12.3g", twall);
            fflush(stdout);

            printf("\n");
        }
    }

    arb_clear(x);
    arb_clear(y);
    arb_clear(res);

    flint_cleanup_master();
    return 0;
}

