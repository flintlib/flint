/* This file is public domain. Author: D.H.J. Polymath. */

#include <string.h>
#include <stdlib.h>
#include <flint/profiler.h>
#include <flint/arb.h>
#include <flint/acb_dirichlet.h>

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
    int automatic;
    slong wp = 0;
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

    /* By default, the library chooses the method (and with the large
     * height method, its working precision) from the height, the count
     * and the target precision: see _acb_dirichlet_hardy_z_zeros_use_platt.
     * An open-ended run counts as many zeros. Don't worry about crossing
     * a threshold, just use the method that is better at the beginning
     * of the run. With -platt, the large height method is used with
     * prec as its working precision (the zeros come out with about
     * prec - log2(t) - 35 bits after the binary point, and at most a
     * height-dependent accuracy: about 187 bits at 1e15); with -noplatt,
     * the Riemann-Siegel method with prec as the target precision.
     */
    automatic = !noplatt && !platt;

    if (prec == -1)
    {
        prec = 64 + fmpz_clog_ui(nstart, 2);
        digits = prec / 3.32192809488736 + 1;
        /* (40 extra bits give the same accuracy as the target precision
           of the other modes) */
        if (platt) prec += 40;
    }

    if (automatic)
    {
        slong len;
        if (fmpz_sgn(requested) < 0 || !fmpz_fits_si(requested))
            len = WORD_MAX;
        else
            len = fmpz_get_si(requested);
        wp = _acb_dirichlet_hardy_z_zeros_use_platt(nstart, len, prec);
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
        else if (noplatt)
        {
            flint_printf("method: noplatt (good for small heights, "
                         "few consecutive zeros, and greater precision; "
                         "interprets prec as a goal precision)\n");
        }
        else if (wp != 0)
        {
            flint_printf("method: automatic, large height method at the "
                         "working precision %wd, then Riemann-Siegel "
                         "refinement where needed (prec is the goal "
                         "precision)\n", wp);
        }
        else
        {
            flint_printf("method: automatic, Riemann-Siegel "
                         "(prec is the goal precision)\n");
        }
    }

    p = _arb_vec_init(buffersize);
    fmpz_set(n, nstart);
    fmpz_zero(count);

    TIMEIT_ONCE_START;

    /* The Riemann-Siegel method or the automatic choice. */
    if (!platt)
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
            /* (the first zeros early with Riemann-Siegel; with the
               large height method, as many as possible at once) */
            if (wp == 0 && fmpz_cmp_si(iter, 30) < 0)
            {
                num = FLINT_MIN(1 << fmpz_get_si(iter), num);
            }
            if (noplatt)
                _acb_dirichlet_hardy_z_zeros_rs(p, n, num, prec);
            else
                acb_dirichlet_hardy_z_zeros(p, n, num, prec);
            print_zeros(p, n, num, digits);
            fmpz_add_si(n, n, num);
            fmpz_add_si(count, count, num);
            fmpz_add_si(iter, iter, 1);
        }
        fmpz_clear(iter);
    }

    /* The large height method at the working precision prec. */
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

    TIMEIT_ONCE_STOP;
    print_memory_usage();

    _arb_vec_clear(p, buffersize);
    fmpz_clear(nstart);
    fmpz_clear(n);
    fmpz_clear(requested);
    fmpz_clear(count);

    flint_cleanup_master();
    return 0;
}
