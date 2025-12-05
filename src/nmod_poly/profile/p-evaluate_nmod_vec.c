/*
   Copyright (C) 2025 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "stdlib.h"  /* for atoi */
#include "nmod.h"
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_vec.h"
#include <string.h>

#define __NB_ITER 10

/*------------------------------------------------------------*/
/* finds an element of order at least n                       */
/* returns 0 if not found                                     */
/*------------------------------------------------------------*/
static long nmod_find_root(slong n, nmod_t mod)
{
    ulong attempts = 0;
    for (ulong q = 2; q < mod.n; q++)
    {
        slong k = 1;
        slong qk = q;
        while (qk != 1 && k < n)
        {
            qk = nmod_mul(qk, q, mod);
            k++;
        }
        if (qk != 1)
        {
            return q;
        }
        attempts += 1;
        if (attempts >= 10)
            return 0;
    }
    return 0;
}

typedef struct
{
    flint_bitcnt_t bits;
    slong length;
    slong npoints;
    slong npoints_precomp;
} info_t;

/* basic iter */
void sample_nmod_vec_iter(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints = info->npoints;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(npoints);
    nn_ptr vals = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(poly, state, length, mod);
        _nmod_vec_rand(pts, state, npoints, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_evaluate_nmod_vec_iter(vals, poly, length, pts, npoints, mod);
        prof_stop();
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* fast, tree */
void sample_nmod_vec_fast(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints = info->npoints;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(npoints);
    nn_ptr vals = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(poly, state, length, mod);
        _nmod_vec_rand(pts, state, npoints, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_evaluate_nmod_vec_fast(vals, poly, length, pts, npoints, mod);
        prof_stop();
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* fast, tree, not counting precomputations */
void sample_nmod_vec_fast_precomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints = info->npoints;
    slong npoints_precomp = info->npoints;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(npoints_precomp);
    nn_ptr vals = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(poly, state, length, mod);
        _nmod_vec_rand(pts, state, npoints_precomp, mod);

        nn_ptr * tree = _nmod_poly_tree_alloc(npoints_precomp);
        _nmod_poly_tree_build(tree, pts, npoints_precomp, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_evaluate_nmod_vec_fast_precomp(vals, poly, length, tree, npoints, mod);
        prof_stop();

        _nmod_poly_tree_free(tree, npoints_precomp);
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* tree, counting only precomputations */
void sample_nmod_vec_fast_onlyprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong npoints_precomp = info->npoints_precomp;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(npoints_precomp);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(pts, state, npoints_precomp, mod);

        nn_ptr * tree = _nmod_poly_tree_alloc(npoints_precomp);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_tree_build(tree, pts, npoints_precomp, mod);
        prof_stop();

        _nmod_poly_tree_free(tree, npoints_precomp);
    }

    _nmod_vec_clear(pts);

    FLINT_TEST_CLEAR(state);
}

/* geometric: basic iter */
void sample_nmod_vec_geom_iter(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints = info->npoints;

    FLINT_TEST_INIT(state);

    nn_ptr vals = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(poly, state, length, mod);
        ulong r = nmod_find_root(2*npoints, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_evaluate_geometric_nmod_vec_iter(vals, poly, length, r, npoints, mod);
        prof_stop();
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, with precomp */
void sample_nmod_vec_geom_fast(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints = info->npoints;

    FLINT_TEST_INIT(state);

    nn_ptr vals = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(poly, state, length, mod);
        ulong r = nmod_find_root(2*npoints, mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_evaluate_geometric_nmod_vec_fast(vals, poly, length, r, npoints, mod);
        prof_stop();
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, not counting precomputations */
void sample_nmod_vec_geom_fast_precomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints = info->npoints;
    slong npoints_precomp = info->npoints_precomp;

    FLINT_TEST_INIT(state);

    nn_ptr vals = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(poly, state, length, mod);
        ulong r = nmod_find_root(2*npoints_precomp, mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        nmod_geometric_progression_t G;
        nmod_geometric_progression_init(G, r, npoints_precomp, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(vals, poly, length, G, npoints, mod);
        prof_stop();

        nmod_geometric_progression_clear(G);
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, counting only precomputations */
void sample_nmod_vec_geom_fast_onlyprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong npoints_precomp = info->npoints_precomp;

    FLINT_TEST_INIT(state);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        ulong r = nmod_find_root(2*npoints_precomp, mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        nmod_geometric_progression_t G;

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_geometric_progression_init(G, r, npoints_precomp, mod);
        prof_stop();

        nmod_geometric_progression_clear(G);
    }

    FLINT_TEST_CLEAR(state);
}

int main(int argc, char * argv[])
{
    if (argc > 1 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        flint_printf("Usage: %s -h for this help, or\n"
                     "       %s [func] [pre_fact] [short]\n"
                     " Optional arguments (if one is provided, previous ones must as well)\n"
                     "       [func] is optional (default 0)\n"
                     "              0 -> all\n"
                     "              1 -> general points only\n"
                     "              2 -> geometric points only\n"
                     "       [pre_fact] is optional (default 1)\n"
                     "              positive integer; precomputation will be on pre_fact * nploints\n"
                     "       [short] is optional (default 0)\n"
                     "               0 -> short bench, 1 -> full bench\n",
                     argv[0], argv[0]);
        return 0;
    }

    const slong func_bench = (argc >= 2) ? atoi(argv[1]) : 0;

    const slong pre_fact = (argc >= 3) ? atoi(argv[2]) : 1;

    flint_printf("[pre_fact provided] -> will build geometric progression with `pre_fact * points` points\n");
    flint_printf("[pre_fact provided] -> note: time for general points w/ tree is for tree with `points` points\n");
    flint_printf("[pre_fact provided] -> note: time for general points precomputation is for tree with `pre_fact * points` points\n");

    /* number of lens differs: long test / short test */
    const slong nb_lens = (argc >= 4 && atoi(argv[3]) != 0) ? 25 : 21;
    slong lengths[] = {1, 2, 3, 4, 6,
                       8, 10, 12, 16, 20,
                       30, 45, 70, 100, 200,
                       400, 800, 1600, 3200, 6400,
                       12800, 25600, 51200, 102400, 204800};

    double tmp;
    double time_iter;
    double time_fast;
    double time_fast_precomp;
    double time_fast_onlyprecomp;
    double time_geom_iter;
    double time_geom_fast;
    double time_geom_fast_precomp;
    double time_geom_fast_onlyprecomp;

    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: measurements in ms\n");

    for (i = 62; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("==== nbits = %ld====\n", i);

        for (int npoints_factor = 1; npoints_factor < 5; npoints_factor += 1)
        {
            if (func_bench == 0)  /* bench all */
            {
                flint_printf("==== nb eval points == %d * (poly length) ====\n", npoints_factor);
                flint_printf("len\tpoints |        GENERAL POINTS         |       GEOMETRIC PROGRESSION     \n");
                flint_printf("len\tpoints |iter\tfast\tw/ tree\ttree   |iter\tfast\tw/ prec\tprecomp\n");
            }
            else if (func_bench == 1)  /* general only */
            {
                flint_printf("==== nb eval points == %d * (poly length) ====\n", npoints_factor);
                flint_printf("len\tpoints |        GENERAL POINTS         |\n");
                flint_printf("len\tpoints |iter\tfast\tw/ tree\ttree   |\n");
            }
            else if (func_bench == 2)  /* geometric only */
            {
                flint_printf("==== nb eval points == %d * (poly length) ====\n", npoints_factor);
                flint_printf("len\tpoints |       GEOMETRIC PROGRESSION     \n");
                flint_printf("len\tpoints |iter\tfast\tw/ prec\tprecomp\n");
            }

            for (int len = 0; len < nb_lens; ++len)
            {
                /* cycles / limb (if CLOCK_SCALE_FACTOR set correctly) */
                /* const double fac = npoints[len] * (double)FLINT_CLOCK_SCALE_FACTOR; */
                /* time in s */
                /* const double fac = 1000000. * __NB_ITER; */
                /* time in ms */
                const double fac = 1. * __NB_ITER;
                info.npoints = npoints_factor * lengths[len];
                info.npoints_precomp = pre_fact * npoints_factor * lengths[len];
                info.length = lengths[len];

                if (func_bench == 0 || func_bench == 1)
                {
                    if (info.npoints <= 1024)
                        prof_repeat(&time_iter, &tmp, sample_nmod_vec_iter, (void *) &info);
                    else
                        time_iter = 0.0;

                    prof_repeat(&time_fast, &tmp, sample_nmod_vec_fast, (void *) &info);

                    prof_repeat(&time_fast_precomp, &tmp, sample_nmod_vec_fast_precomp, (void *) &info);

                    prof_repeat(&time_fast_onlyprecomp, &tmp, sample_nmod_vec_fast_onlyprecomp, (void *) &info);
                }

                if (func_bench == 0 || func_bench == 2)
                {
                    if (info.npoints <= 1024)
                        prof_repeat(&time_geom_iter, &tmp, sample_nmod_vec_geom_iter, (void *) &info);
                    else
                        time_geom_iter = 0.0;

                    prof_repeat(&time_geom_fast, &tmp, sample_nmod_vec_geom_fast, (void *) &info);

                    prof_repeat(&time_geom_fast_precomp, &tmp, sample_nmod_vec_geom_fast_precomp, (void *) &info);

                    prof_repeat(&time_geom_fast_onlyprecomp, &tmp, sample_nmod_vec_geom_fast_onlyprecomp, (void *) &info);
                }

                if (func_bench == 0)
                {
                    flint_printf("%ld\t%7ld|%.1e\t%.1e\t%.1e\t%.1e|%.1e\t%.1e\t%.1e\t%.1e\n",
                                 info.length, info.npoints,
                                 time_iter/fac, time_fast/fac, time_fast_precomp/fac, time_fast_onlyprecomp/fac,
                                 time_geom_iter/fac, time_geom_fast/fac, time_geom_fast_precomp/fac, time_geom_fast_onlyprecomp/fac);
                }
                else if (func_bench == 1)
                {
                    flint_printf("%ld\t%7ld|%.1e\t%.1e\t%.1e\t%.1e\n",
                                 info.length, info.npoints,
                                 time_iter/fac, time_fast/fac, time_fast_precomp/fac, time_fast_onlyprecomp/fac);
                }
                else if (func_bench == 2)
                {
                    flint_printf("%ld\t%7ld|%.1e\t%.1e\t%.1e\t%.1e\n",
                                 info.length, info.npoints,
                                 time_geom_iter/fac, time_geom_fast/fac, time_geom_fast_precomp/fac, time_geom_fast_onlyprecomp/fac);
                }
            }
        }

        flint_printf("\n");
    }

    return 0;
}

#undef __NB_ITER
