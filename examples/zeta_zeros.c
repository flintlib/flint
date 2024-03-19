/* This file is public domain. Author: D.H.J. Polymath. */

#include <string.h>
#include <stdlib.h>
#include <flint/acb_dirichlet.h>
#include <flint/profiler.h>

void print_zeros(arb_srcptr p, const fmpz_t n, slong len, slong digits)
{
    slong i;
    fmpz_t k;
    fmpz_init_set(k, n);
    for (i = 0; i < len; i++)
    {
        fmpz_print(k);
        flint_printf("\t");
        arb_printn(p+i, digits, ARB_STR_NO_RADIUS);
        flint_printf("\n");
        fmpz_add_ui(k, k, 1);
    }
    fmpz_clear(k);
}

void print_help(void)
{
    flint_printf("zeta_zeros [-n n] [-count n] [-prec n] [-digits n] [-threads n] "
                 "[-platt] [-noplatt] [-v] [-verbose] [-h] [-help]\n\n");
    flint_printf("Reports the imaginary parts of consecutive nontrivial zeros "
                 "of the Riemann zeta function.\n");
}

void requires_value(int argc, char *argv[], slong i)
{
    if (i == argc-1)
    {
        flint_printf("the argument %s requires a value\n", argv[i]);
        flint_abort();
    }
}

void invalid_value(char *argv[], slong i)
{
    flint_printf("invalid value for the argument %s: %s\n", argv[i], argv[i+1]);
    flint_abort();
}

int main(int argc, char *argv[])
{
    const slong max_buffersize = 30000;
    int verbose = 0;
    int platt = 0;
    int noplatt = 0;
    slong i, buffersize, prec, digits;
    fmpz_t requested, count, nstart, n;
    arb_ptr p;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help"))
        {
            print_help();
            return 0;
        }
    }

    fmpz_init(requested);
    fmpz_init(count);
    fmpz_init(nstart);
    fmpz_init(n);

    fmpz_one(nstart);
    fmpz_set_si(requested, -1);
    buffersize = max_buffersize;
    prec = -1;
    digits = 2;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-noplatt"))
        {
            noplatt = 1;
        }
        else if (!strcmp(argv[i], "-platt"))
        {
            platt = 1;
        }
        else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "-verbose"))
        {
            verbose = 1;
        }
        else if (!strcmp(argv[i], "-threads"))
        {
            slong threads;
            requires_value(argc, argv, i);
            threads = atol(argv[i+1]);
            if (threads < 1)
            {
                invalid_value(argv, i);
            }
            flint_set_num_threads(threads);
        }
        else if (!strcmp(argv[i], "-n"))
        {
            requires_value(argc, argv, i);
            if (fmpz_set_str(nstart, argv[i+1], 10) || fmpz_sgn(nstart) < 1)
            {
                invalid_value(argv, i);
            }
        }
        else if (!strcmp(argv[i], "-count"))
        {
            requires_value(argc, argv, i);
            if (fmpz_set_str(requested, argv[i+1], 10) ||
                fmpz_sgn(requested) < 1)
            {
                invalid_value(argv, i);
            }
            if (fmpz_cmp_si(requested, buffersize) < 0)
            {
                buffersize = fmpz_get_si(requested);
            }
        }
        else if (!strcmp(argv[i], "-prec"))
        {
            requires_value(argc, argv, i);
            prec = atol(argv[i+1]);
            digits = prec / 3.32192809488736 + 1;
            if (prec < 2)
            {
                invalid_value(argv, i);
            }
        }
        else if (!strcmp(argv[i], "-digits"))
        {
            requires_value(argc, argv, i);
            digits = atol(argv[i+1]);
            prec = digits * 3.32192809488736 + 3;
            if (prec < 2)
            {
                invalid_value(argv, i);
            }
        }
    }

    if (platt && noplatt)
    {
        flint_printf("conflicting arguments platt and noplatt\n");
        flint_abort();
    }

    if (platt && fmpz_cmp_si(nstart, 10000) < 0)
    {
        flint_printf("this implementation of the platt algorithm "
                     "is not valid below the 10000th zero\n");
        flint_abort();
    }

    /* Above n~1e15 the large height method is better.
     * Above n~1e11 the large height method is better for more than 100 zeros.
     * Don't worry about crossing the threshold, just use the method
     * that is better at the beginning of the run.
     */
    if (!noplatt && !platt)
    {
        fmpz_t threshold;
        fmpz_init(threshold);
        fmpz_set_si(threshold, 10);
        fmpz_pow_ui(threshold, threshold, 15);
        if (fmpz_sgn(requested) < 0 || fmpz_cmp_si(requested, 100) > 0)
        {
            fmpz_set_si(threshold, 10);
            fmpz_pow_ui(threshold, threshold, 11);
        }
        if (fmpz_cmp(nstart, threshold) < 0)
        {
            noplatt = 1;
        }
        else
        {
            platt = 1;
        }
        fmpz_clear(threshold);
    }

    if (prec == -1)
    {
        prec = 64 + fmpz_clog_ui(nstart, 2);
        if (platt) prec *= 2;
        digits = prec / 3.32192809488736 + 1;
    }

    if (verbose)
    {
        flint_printf("n: "); fmpz_print(nstart); flint_printf("\n");
        flint_printf("count: "); fmpz_print(requested); flint_printf("\n");
        flint_printf("threads: %wd\n", flint_get_num_threads());
        flint_printf("prec: %wd\n", prec);
        if (platt)
        {
            flint_printf("method: platt (good for large heights, "
                         "many consecutive zeros, and lower precision; "
                         "interprets prec as a working precision)\n");
        }
        else
        {
            flint_printf("method: noplatt (good for small heights, "
                         "few consecutive zeros, and greater precision; "
                         "interprets prec as a goal precision)\n");
        }
    }

    p = _arb_vec_init(buffersize);
    fmpz_set(n, nstart);
    fmpz_zero(count);

    TIMEIT_ONCE_START

    /* This is the low height method. */
    if (noplatt)
    {
        fmpz_t iter;
        fmpz_init(iter);
        while (fmpz_sgn(requested) < 0 || fmpz_cmp(count, requested) < 0)
        {
            slong num = buffersize;
            if (fmpz_sgn(requested) >= 0)
            {
                fmpz_t remaining;
                fmpz_init(remaining);
                fmpz_sub(remaining, requested, count);
                if (fmpz_cmp_si(remaining, num) < 0)
                {
                    num = fmpz_get_si(remaining);
                }
                fmpz_clear(remaining);
            }
            if (fmpz_cmp_si(iter, 30) < 0)
            {
                num = FLINT_MIN(1 << fmpz_get_si(iter), num);
            }
            acb_dirichlet_hardy_z_zeros(p, n, num, prec);
            print_zeros(p, n, num, digits);
            fmpz_add_si(n, n, num);
            fmpz_add_si(count, count, num);
            fmpz_add_si(iter, iter, 1);
        }
        fmpz_clear(iter);
    }

    /* This is the large height method. */
    if (platt)
    {
        while (fmpz_sgn(requested) < 0 || fmpz_cmp(count, requested) < 0)
        {
            slong found;
            slong num = buffersize;
            if (fmpz_sgn(requested) >= 0)
            {
                fmpz_t remaining;
                fmpz_init(remaining);
                fmpz_sub(remaining, requested, count);
                if (fmpz_cmp_si(remaining, num) < 0)
                {
                    num = fmpz_get_si(remaining);
                }
                fmpz_clear(remaining);
            }
            found = acb_dirichlet_platt_local_hardy_z_zeros(p, n, num, prec);
            if (!found)
            {
                flint_printf("Failed to find some zeros.\n");
                flint_printf("Maybe prec is not high enough or something "
                             "is wrong with the internal tuning parameters.");
                flint_abort();
            }
            print_zeros(p, n, found, digits);
            fmpz_add_si(n, n, found);
            fmpz_add_si(count, count, found);
        }
    }

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    _arb_vec_clear(p, buffersize);
    fmpz_clear(nstart);
    fmpz_clear(n);
    fmpz_clear(requested);
    fmpz_clear(count);

    flint_cleanup_master();
    return 0;
}
