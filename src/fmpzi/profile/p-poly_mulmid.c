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
    Candidate algorithms for a future _fmpzi_poly_mulmid dispatcher.

    Usage: p-poly_mulmid [options]

      -s SHAPE   mul | sqr | mullow | mulhigh | mid | unbal | all
                 (as in fmpz_poly/profile/p-mulmid_SS; default mul)
      -i IMAG    balanced  real and imaginary parts of the same size
                           (default)
                 quarter   imaginary parts of the first operand have
                           a quarter of the bits
                 real      the first operand is real
      -t T       number of threads (default 1)
      -b B       largest size of the real and imaginary parts in bits
                 (default 262144)
      -n N       largest length n (default 4096)
      -w W       skip cells with n * bits > W (default 2^24)
      -g G       lengths per octave (default 2)
      -q Q       skip the quadratic algorithms when n > Q (default 128)
      -v         also print the timings of every cell (to stderr)

    The candidates are

      C  _gr_poly_mulmid_classical over fmpzi
      F  _fmpzi_poly_mulmid_classical_fft_small (when it accepts the input)
      4  four real _fmpz_poly_mulmid products (complex reorder, karatsuba = 0)
      3  three real products (complex reorder, karatsuba = 1)
      P  _fmpzi_poly_mulmid_SS_pair
      E  _fmpzi_poly_mulmid_SS_single (single embedding)
      K  _fmpzi_poly_mulmid_KS

    (_fmpzi_poly_mulmid_SS chooses between P and E).

    and the reference is the dispatcher _fmpzi_poly_mulmid (which the
    gr fmpzi ring uses). For each shape it prints two matrices with the
    sizes as rows and lengths as columns: the fastest candidate, and
    t(fastest) / t(dispatcher), i.e. what an ideal choice among the
    candidates would gain over the dispatcher -- values well below 1
    mark cells where its boundaries are off (it always runs one of the
    candidates, so values above 1 are timing noise).
    A candidate is not
    repeated once a single call exceeds 8 times the best time so far.
    All results are checked against each other.
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpzi.h"
#include "gr.h"
#include "gr_poly.h"
#include "thread_support.h"
#include "profiler.h"

#define NALG 7

static const char alg_char[NALG] = { 'C', 'F', '4', '3', 'P', 'E', 'K' };

static const char * shapes[] = { "mul", "sqr", "mullow", "mulhigh", "mid", "unbal" };

typedef struct
{
    int alg;                /* 0..NALG-1, or -1 for the dispatcher */
    fmpzi_struct * res;
    const fmpzi_struct * a;
    const fmpzi_struct * b;
    slong len1, len2, nlo, nhi;
    gr_ctx_struct * ctx;
    gr_ctx_struct * zctx;
    int ok;                 /* the candidate accepted the input */
}
call_struct;

static void
do_call(call_struct * x)
{
    int status = GR_SUCCESS;

    switch (x->alg)
    {
        case 0:
            status = _gr_poly_mulmid_classical(x->res, x->a, x->len1, x->b, x->len2,
                                               x->nlo, x->nhi, x->ctx);
            break;
        case 1:
            x->ok = _fmpzi_poly_mulmid_classical_fft_small(x->res, x->a, x->len1,
                                               x->b, x->len2, x->nlo, x->nhi);
            break;
        case 2:
        case 3:
            status = _gr_poly_mulmid_complex_reorder(x->res, x->a, x->len1, x->b, x->len2,
                                     x->nlo, x->nhi, x->alg == 3, x->ctx, x->zctx);
            break;
        case 4:
            _fmpzi_poly_mulmid_SS_pair(x->res, x->a, x->len1, x->b, x->len2, x->nlo, x->nhi);
            break;
        case 5:
            _fmpzi_poly_mulmid_SS_single(x->res, x->a, x->len1, x->b, x->len2, x->nlo, x->nhi);
            break;
        case 6:
            _fmpzi_poly_mulmid_KS(x->res, x->a, x->len1, x->b, x->len2, x->nlo, x->nhi);
            break;
        default:
            _fmpzi_poly_mulmid(x->res, x->a, x->len1, x->b, x->len2, x->nlo, x->nhi);
    }

    if (status != GR_SUCCESS)
        flint_abort();
}

/* wall time in microseconds per call: best of three batches of at least
   20 ms, or a single warm call when that exceeds both cap and 1 ms */
static double
time_call(call_struct * x, double cap)
{
    double best = 1e300;
    slong batch, reps, i;
    timeit_t t;

    x->ok = 1;
    timeit_start_us(t);
    do_call(x);
    timeit_stop_us(t);

    if (!x->ok)
        return -1.0;

    if (t->wall > 200000)
        return (double) t->wall;

    /* the first call can be inflated by cold caches and lazy
       initialization; decide on the cap from a warm call, and only for
       calls long enough that repeating them would cost real time */
    timeit_start_us(t);
    do_call(x);
    timeit_stop_us(t);

    if (t->wall > 1000 && (double) t->wall > cap)
        return (double) t->wall;

    for (batch = 0; batch < 3; batch++)
    {
        reps = 1;
        for (;;)
        {
            timeit_start_us(t);
            for (i = 0; i < reps; i++)
                do_call(x);
            timeit_stop_us(t);

            if (t->wall >= 20000)
                break;
            reps *= (t->wall < 1000) ? 10 : 2;
        }
        best = FLINT_MIN(best, (double) t->wall / reps);
    }

    return best;
}

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

static fmpzi_struct *
vec_init(slong len)
{
    slong i;
    fmpzi_struct * v = flint_malloc(len * sizeof(fmpzi_struct));
    for (i = 0; i < len; i++)
        fmpzi_init(v + i);
    return v;
}

static void
vec_clear(fmpzi_struct * v, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        fmpzi_clear(v + i);
    flint_free(v);
}

static void
profile_shape(int shape, int imag, slong maxbits, slong maxlen, double maxwork,
              int grid, slong maxquad, int verbose, flint_rand_t state)
{
    gr_ctx_t ctx, zctx;
    slong lens[64], bitsv[64], nlens = 0, nbits = 0;
    slong li, bi, n, bits, k;
    char * winner;
    double * gain;

    gr_ctx_init_fmpzi(ctx);
    gr_ctx_init_fmpz(zctx);

    for (k = 0; nlens < 64; k++)
    {
        n = (slong) (4.0 * pow(2.0, (double) k / grid) + 0.5);
        if (n > maxlen)
            break;
        if (nlens == 0 || n != lens[nlens - 1])
            lens[nlens++] = n;
    }

    for (bits = 64; bits <= maxbits && nbits < 64; bits *= 2)
        bitsv[nbits++] = bits;

    winner = flint_calloc(nbits * nlens, 1);
    gain = flint_calloc(nbits * nlens, sizeof(double));

    flint_printf("\nshape %s, imaginary parts %s\n", shapes[shape],
        imag == 0 ? "balanced" : imag == 1 ? "quarter size (first operand)" : "zero (first operand)");

    for (bi = 0; bi < nbits; bi++)
    {
        bits = bitsv[bi];

        for (li = 0; li < nlens; li++)
        {
            slong len1, len2, nlo, nhi, i, alg;
            fmpzi_struct * a, * b, * ref, * r;
            call_struct x;
            double t[NALG], tgr, best = 1e300;
            int bestalg = -1;

            n = lens[li];

            if ((double) n * (double) bits > maxwork)
                continue;

            shape_params(shape, n, &len1, &len2, &nlo, &nhi);

            a = vec_init(len1);
            b = (shape == 1) ? a : vec_init(len2);
            ref = vec_init(nhi - nlo);
            r = vec_init(nhi - nlo);

            for (i = 0; i < len1; i++)
            {
                fmpz_randbits(fmpzi_realref(a + i), state, bits);
                if (imag == 0)
                    fmpz_randbits(fmpzi_imagref(a + i), state, bits);
                else if (imag == 1)
                    fmpz_randbits(fmpzi_imagref(a + i), state, FLINT_MAX(bits / 4, 1));
            }
            if (shape != 1)
            {
                for (i = 0; i < len2; i++)
                {
                    fmpz_randbits(fmpzi_realref(b + i), state, bits);
                    fmpz_randbits(fmpzi_imagref(b + i), state, bits);
                }
            }

            x.a = a; x.b = b; x.len1 = len1; x.len2 = len2;
            x.nlo = nlo; x.nhi = nhi; x.ctx = ctx; x.zctx = zctx;

            /* reference: the dispatcher */
            x.alg = -1;
            x.res = ref;
            tgr = time_call(&x, 1e300);

            /* fast candidates first, so that the cap is useful */
            for (k = NALG - 1; k >= 0; k--)
            {
                alg = k;
                t[alg] = -1.0;

                if ((alg == 0 || alg == 1) && FLINT_MIN(len1, len2) > maxquad)
                    continue;

                x.alg = alg;
                x.res = r;
                t[alg] = time_call(&x, 8 * FLINT_MIN(best, tgr));

                if (t[alg] < 0)
                    continue;

                for (i = 0; i < nhi - nlo; i++)
                {
                    if (!fmpzi_equal(r + i, ref + i))
                    {
                        flint_printf("\nFAIL: %c differs (shape %s, bits %wd, n %wd)\n",
                                     alg_char[alg], shapes[shape], bits, n);
                        flint_abort();
                    }
                }

                if (t[alg] < best)
                {
                    best = t[alg];
                    bestalg = alg;
                }
            }

            winner[bi * nlens + li] = alg_char[bestalg];
            gain[bi * nlens + li] = best / tgr;

            if (verbose)
            {
                flint_fprintf(stderr, "%s bits=%wd n=%wd  D %.4g", shapes[shape], bits, n, tgr);
                for (alg = 0; alg < NALG; alg++)
                {
                    if (t[alg] < 0)
                        flint_fprintf(stderr, "  %c -", alg_char[alg]);
                    else
                        flint_fprintf(stderr, "  %c %.4g", alg_char[alg], t[alg]);
                }
                flint_fprintf(stderr, " us\n");
            }

            vec_clear(a, len1);
            if (shape != 1)
                vec_clear(b, len2);
            vec_clear(ref, nhi - nlo);
            vec_clear(r, nhi - nlo);
        }
    }

    flint_printf("fastest candidate (C classical, F classical fft_small, 4/3 real products, P SS pair, E SS single, K KS)\n");
    flint_printf("%8s", "bits\\n");
    for (li = 0; li < nlens; li++)
        flint_printf(" %5wd", lens[li]);
    flint_printf("\n");
    for (bi = 0; bi < nbits; bi++)
    {
        flint_printf("%8wd", bitsv[bi]);
        for (li = 0; li < nlens; li++)
            flint_printf("     %c", winner[bi * nlens + li] ? winner[bi * nlens + li] : ' ');
        flint_printf("\n");
    }

    flint_printf("t(fastest) / t(_fmpzi_poly_mulmid)\n");
    flint_printf("%8s", "bits\\n");
    for (li = 0; li < nlens; li++)
        flint_printf(" %5wd", lens[li]);
    flint_printf("\n");
    for (bi = 0; bi < nbits; bi++)
    {
        flint_printf("%8wd", bitsv[bi]);
        for (li = 0; li < nlens; li++)
        {
            if (winner[bi * nlens + li])
                flint_printf(" %5.2f", gain[bi * nlens + li]);
            else
                flint_printf(" %5s", "");
        }
        flint_printf("\n");
    }

    fflush(stdout);
    flint_free(winner);
    flint_free(gain);
    gr_ctx_clear(ctx);
    gr_ctx_clear(zctx);
}

int main(int argc, char * argv[])
{
    slong maxbits = WORD(1) << 18, maxlen = 4096, threads = 1, maxquad = 128;
    double maxwork = 16777216.0;
    int shape = 0, imag = 0, verbose = 0, all = 0, grid = 2, i;
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
                for (k = 0; k < 6; k++)
                    if (!strcmp(s, shapes[k]))
                        shape = k;
        }
        else if (!strcmp(argv[i], "-i") && i + 1 < argc)
        {
            const char * s = argv[++i];
            imag = !strcmp(s, "quarter") ? 1 : !strcmp(s, "real") ? 2 : 0;
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
        else if (!strcmp(argv[i], "-q") && i + 1 < argc)
            maxquad = atol(argv[++i]);
        else if (!strcmp(argv[i], "-v"))
            verbose = 1;
        else
        {
            flint_printf("usage: %s [-s mul|sqr|mullow|mulhigh|mid|unbal|all]"
                         " [-i balanced|quarter|real] [-t threads] [-b maxbits]"
                         " [-n maxlen] [-w maxwork] [-g per_octave] [-q maxquad] [-v]\n",
                         argv[0]);
            return 1;
        }
    }

    flint_set_num_threads(threads);
    flint_printf("threads = %wd\n", threads);

    if (all)
    {
        for (shape = 0; shape < 6; shape++)
            profile_shape(shape, imag, maxbits, maxlen, maxwork, grid, maxquad,
                          verbose, state);
    }
    else
        profile_shape(shape, imag, maxbits, maxlen, maxwork, grid, maxquad,
                      verbose, state);

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
