/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Where does _fmpz_poly_mulmid_SS beat the dispatcher _fmpz_poly_mulmid?

    Usage: p-mulmid_SS [options]

      -s SHAPE   mul     n x n -> [0, 2n-1)                 (default)
                 sqr     squaring, n -> [0, 2n-1)
                 mullow  n x n -> [0, n)
                 mulhigh n x n -> [n-1, 2n-1)
                 mid     (2n-1) x n -> [n-1, 2n-1)   (Newton middle product)
                 unbal   4n x n -> [0, 5n-1)
                 all     each of the above in turn
      -t T       number of threads (default 1)
      -b B       largest coefficient size in bits (default 1048576)
      -n N       largest length n (default 8192)
      -w W       skip cells with n * bits > W (default 2^26)
      -g G       lengths per octave (default 2: 4, 6, 8, 11, 16, ...)
      -v         also print one line of timings per cell (to stderr)

    For each shape, prints a matrix of t(SS) / t(_fmpz_poly_mulmid) with
    coefficient sizes as rows and lengths as columns; '*' marks cells
    where SS is at least 5% faster, '=' cells within 2% (typically
    where the dispatcher already chose SS). The lengths
    include non-powers of two, since SS (and its cyclic middle-product
    transform) depends on where 2n - 1 falls between powers of two.
    The coefficients are uniformly random signed values of the given
    size; the results of SS and the dispatcher are compared once per
    cell.
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "thread_support.h"
#include "profiler.h"

typedef void (*mulmid_fn)(fmpz *, const fmpz *, slong, const fmpz *, slong, slong, slong);

typedef struct
{
    mulmid_fn f;
    fmpz * res;
    const fmpz * a;
    const fmpz * b;
    slong len1, len2, nlo, nhi;
}
mulmid_arg_struct;

/* best of three batches, each repeated to at least 20 ms; wall time in
   microseconds per call (wall time, since the transforms may use
   several threads) */
static double
time_call(mulmid_arg_struct * x)
{
    double best = 1e300;
    slong batch, reps, i;
    timeit_t t;

    timeit_start_us(t);
    x->f(x->res, x->a, x->len1, x->b, x->len2, x->nlo, x->nhi);
    timeit_stop_us(t);

    /* a single slow call is accurate enough */
    if (t->wall > 200000)
        return (double) t->wall;

    for (batch = 0; batch < 3; batch++)
    {
        reps = 1;
        for (;;)
        {
            timeit_start_us(t);
            for (i = 0; i < reps; i++)
                x->f(x->res, x->a, x->len1, x->b, x->len2, x->nlo, x->nhi);
            timeit_stop_us(t);

            if (t->wall >= 20000)
                break;
            reps *= (t->wall < 1000) ? 10 : 2;
        }
        best = FLINT_MIN(best, (double) t->wall / reps);
    }

    return best;
}

static const char * shapes[] = { "mul", "sqr", "mullow", "mulhigh", "mid", "unbal" };

static void
shape_params(int shape, slong n, slong * len1, slong * len2, slong * nlo, slong * nhi)
{
    switch (shape)
    {
        case 0: case 1:
            *len1 = n; *len2 = n; *nlo = 0; *nhi = 2 * n - 1; break;
        case 2:
            *len1 = n; *len2 = n; *nlo = 0; *nhi = n; break;
        case 3:
            *len1 = n; *len2 = n; *nlo = n - 1; *nhi = 2 * n - 1; break;
        case 4:
            *len1 = 2 * n - 1; *len2 = n; *nlo = n - 1; *nhi = 2 * n - 1; break;
        default:
            *len1 = 4 * n; *len2 = n; *nlo = 0; *nhi = 5 * n - 1; break;
    }
}

static void
profile_shape(int shape, slong maxbits, slong maxlen, double maxwork, int grid,
              int verbose, flint_rand_t state)
{
    slong lens[64], bitsv[64], nlens = 0, nbits = 0;
    slong li, bi, n, bits;

    /* lengths from 4, geometrically with 'grid' points per octave */
    {
        slong k;
        for (k = 0; nlens < 64; k++)
        {
            n = (slong) (4.0 * pow(2.0, (double) k / grid) + 0.5);
            if (n > maxlen)
                break;
            if (nlens == 0 || n != lens[nlens - 1])
                lens[nlens++] = n;
        }
    }

    for (bits = 64; bits <= maxbits; bits *= 2)
        bitsv[nbits++] = bits;

    flint_printf("\nshape %s: t(SS) / t(_fmpz_poly_mulmid), rows = bits, columns = n\n",
                 shapes[shape]);
    flint_printf("%8s", "bits\\n");
    for (li = 0; li < nlens; li++)
        flint_printf(" %6wd", lens[li]);
    flint_printf("\n");

    for (bi = 0; bi < nbits; bi++)
    {
        bits = bitsv[bi];
        flint_printf("%8wd", bits);

        for (li = 0; li < nlens; li++)
        {
            slong len1, len2, nlo, nhi, i;
            fmpz * a, * b, * r1, * r2;
            mulmid_arg_struct x;
            double tdef, tss, ratio;

            n = lens[li];

            if ((double) n * (double) bits > maxwork)
            {
                flint_printf(" %6s", "");
                continue;
            }

            shape_params(shape, n, &len1, &len2, &nlo, &nhi);

            a = _fmpz_vec_init(len1);
            b = (shape == 1) ? a : _fmpz_vec_init(len2);
            r1 = _fmpz_vec_init(nhi - nlo);
            r2 = _fmpz_vec_init(nhi - nlo);

            for (i = 0; i < len1; i++)
                fmpz_randbits(a + i, state, bits);
            if (shape != 1)
                for (i = 0; i < len2; i++)
                    fmpz_randbits(b + i, state, bits);

            x.a = a; x.b = b; x.len1 = len1; x.len2 = len2;
            x.nlo = nlo; x.nhi = nhi;

            x.f = _fmpz_poly_mulmid;
            x.res = r1;
            tdef = time_call(&x);

            x.f = _fmpz_poly_mulmid_SS;
            x.res = r2;
            tss = time_call(&x);

            if (!_fmpz_vec_equal(r1, r2, nhi - nlo))
            {
                flint_printf("\nFAIL: results differ (shape %s, bits %wd, n %wd)\n",
                             shapes[shape], bits, n);
                flint_abort();
            }

            ratio = tss / tdef;
            flint_printf(" %5.2f%c", ratio,
                ratio < 0.95 ? '*' : (ratio < 1.02 && ratio > 0.98) ? '=' : ' ');

            if (verbose)
                flint_fprintf(stderr, "%s bits=%wd n=%wd len1=%wd len2=%wd window=[%wd,%wd)"
                             "  default %.4g us  SS %.4g us\n",
                             shapes[shape], bits, n, len1, len2, nlo, nhi, tdef, tss);

            fflush(stdout);

            _fmpz_vec_clear(a, len1);
            if (shape != 1)
                _fmpz_vec_clear(b, len2);
            _fmpz_vec_clear(r1, nhi - nlo);
            _fmpz_vec_clear(r2, nhi - nlo);
        }

        flint_printf("\n");
    }
}

int main(int argc, char * argv[])
{
    slong maxbits = WORD(1) << 20, maxlen = 8192, threads = 1;
    double maxwork = 67108864.0;
    int shape = 0, verbose = 0, all = 0, grid = 2, i;
    FLINT_TEST_INIT(state);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-s") && i + 1 < argc)
        {
            const char * s = argv[++i];
            int k;
            if (!strcmp(s, "all"))
                all = 1;
            else
            {
                for (k = 0; k < 6; k++)
                    if (!strcmp(s, shapes[k]))
                        shape = k;
            }
        }
        else if (!strcmp(argv[i], "-t") && i + 1 < argc)
            threads = atol(argv[++i]);
        else if (!strcmp(argv[i], "-b") && i + 1 < argc)
            maxbits = atol(argv[++i]);
        else if (!strcmp(argv[i], "-n") && i + 1 < argc)
            maxlen = atol(argv[++i]);
        else if (!strcmp(argv[i], "-w") && i + 1 < argc)
            maxwork = atof(argv[++i]);
        else if (!strcmp(argv[i], "-g") && i + 1 < argc)
        {
            grid = atoi(argv[++i]);
            grid = FLINT_MAX(grid, 1);
        }
        else if (!strcmp(argv[i], "-v"))
            verbose = 1;
        else
        {
            flint_printf("usage: %s [-s mul|sqr|mullow|mulhigh|mid|unbal|all] [-t threads]"
                         " [-b maxbits] [-n maxlen] [-w maxwork] [-g per_octave] [-v]\n", argv[0]);
            return 1;
        }
    }

    flint_set_num_threads(threads);
    flint_printf("threads = %wd\n", threads);

    if (all)
    {
        for (shape = 0; shape < 6; shape++)
            profile_shape(shape, maxbits, maxlen, maxwork, grid, verbose, state);
    }
    else
        profile_shape(shape, maxbits, maxlen, maxwork, grid, verbose, state);

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
