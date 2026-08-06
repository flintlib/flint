/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/profiler.h>

/* There are two versions of the Brent-Salamin AGM iteration, producing
   lower and upper bounds for pi:

(lower)
    a = 1
    b = sqrt(2) / 2
    s = 1
    for k=1..n
        a' = (a+b) / 2
        b' = sqrt(a*b)
        s' = s - 2^(k-1) * (a-b)^2
        a, b, s = a', b', s'
    res = (a+b)^2 / s

(upper)
    a = 1
    b = sqrt(2) / 2
    s = 2*b - 1/2
    for k=1..n
        a' = (a+b) / 2
        b' = sqrt(a*b)
        s' = s - 2^k * (a'-b')^2
        a, b, s = a', b', s'
    res = (a+b)^2 / s

n   (lower)           (upper)
1   3.1405...         3.141680...
2   3.141592646...    3.14159265389...

We implement the upper bound which is very slightly more accurate and admits
the following convenient error bound from Salamin's paper:

    err         <= 2^(n+8) exp(-pi 2^(n+1))
    2^err_2exp  <= (n + 8) - (pi / log(2)) * 2^(n+1)
    pi/log(2)   = 4.532360141...

References:

* Brent - https://maths-people.anu.edu.au/~brent/pd/rpb028.pdf
* Salamin - https://doi.org/10.2307/2005327

*/

/* Using a generous lower bound for pi/log(2), this expression in
   machine arithmetic is sure to give a correct upper bound. */
static slong pi_agm_err_2exp(slong n)
{
    return (n + 8) - 4.53236 * ((double) (UWORD(1) << (n + 1))) + 1;
}

void
arb_pi_agm(arb_t res, slong prec, int verbose)
{
    slong k, n, err_2exp, wp;
    arb_t a, b, s, t;

    arb_init(a);
    arb_init(b);
    arb_init(s);
    arb_init(t);

    /* Heuristically increase working precision to ensure near-correct
       rounding to prec bits. */
    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);

    for (n = 1; ; n++)
    {
        err_2exp = pi_agm_err_2exp(n);
        if (err_2exp <= -wp)
            break;
    }

    arb_one(a);
    arb_sqrt_ui(b, 2, wp);
    arb_mul_2exp_si(b, b, -1);
    arb_set_d(s, -0.25);
    arb_add(s, s, b, wp);
    arb_mul_2exp_si(s, s, 1);

    for (k = 1; k <= n; k++)
    {
        if (verbose)
        {
            flint_printf("Iteration %wd / %wd: error bound ", k, n);
            arb_one(t);
            arb_mul_2exp_si(t, t, pi_agm_err_2exp(k));
            arb_printn(t, 3, ARB_STR_NO_RADIUS);
            flint_printf("\n");
        }

        arb_add(t, a, b, wp);
        arb_mul_2exp_si(t, t, -1);
        arb_mul(b, a, b, wp);
        arb_sqrt(b, b, wp);
        arb_swap(a, t);

        arb_sub(t, a, b, wp);
        arb_sqr(t, t, wp);
        arb_mul_2exp_si(t, t, k);
        arb_sub(s, s, t, wp);
    }

    arb_clear(t);
    arb_add(a, a, b, wp);
    arb_clear(b);
    arb_sqr(a, a, wp);
    arb_div(res, a, s, prec);

    arb_add_error_2exp_si(res, err_2exp);

    arb_clear(a);
    arb_clear(s);
}



int main(int argc, char *argv[])
{
    arb_t x;
    slong i, prec, digits, condense, num_threads;

    if (argc < 2)
    {
        flint_printf("usage: build/examples/pi_agm digits [-condense n] [-threads n]\n");
        flint_printf("compute digits of pi using AGM iteration\n");
        return 1;
    }

    digits = atol(argv[1]);

    num_threads = 1;
    condense = 20;

    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-condense"))
            condense = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
    }

    flint_set_num_threads(num_threads);

    arb_init(x);

    prec = digits * 3.3219280948873623 + 5;

    flint_printf("precision = %wd bits...\n", prec);
    fflush(stdout);

    TIMEIT_ONCE_START;
    arb_pi_agm(x, prec, 1);
    TIMEIT_ONCE_STOP;

    print_memory_usage();

    arb_printn(x, digits, ARB_STR_CONDENSE * condense);

    flint_printf("\n");

    arb_clear(x);
    flint_cleanup_master();
    return 0;
}

