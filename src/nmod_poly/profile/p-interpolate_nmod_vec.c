/*
   Copyright (C) 2026 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  /* for atoi */
#include <string.h>

#include "nmod.h"
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_vec.h"

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
    slong npoints_precomp;
} info_t;

/* Newton */
void sample_nmod_vec_newton(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(length);
    nn_ptr vals = _nmod_vec_init(length);
    nn_ptr poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(vals, state, length, mod);
        _nmod_vec_rand(pts, state, length, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_interpolate_nmod_vec_newton(poly, pts, vals, length, mod);
        prof_stop();
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);
    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

/* barycentric */
void sample_nmod_vec_barycentric(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(length);
    nn_ptr vals = _nmod_vec_init(length);
    nn_ptr poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(vals, state, length, mod);
        _nmod_vec_rand(pts, state, length, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_interpolate_nmod_vec_barycentric(poly, pts, vals, length, mod);
        prof_stop();
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);
    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

/* fast */
void sample_nmod_vec_fast(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(length);
    nn_ptr vals = _nmod_vec_init(length);
    nn_ptr poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(vals, state, length, mod);
        _nmod_vec_rand(pts, state, length, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_interpolate_nmod_vec_fast(poly, pts, vals, length, mod);
        prof_stop();
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);
    _nmod_vec_clear(poly);

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

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(length);
    nn_ptr vals = _nmod_vec_init(length);
    nn_ptr poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(vals, state, length, mod);
        _nmod_vec_rand(pts, state, length, mod);

        nn_ptr * tree = _nmod_poly_tree_alloc(length);
        _nmod_poly_tree_build(tree, pts, length, mod);

        nn_ptr w = _nmod_vec_init(length);
        _nmod_poly_interpolation_weights(w, tree, length, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_interpolate_nmod_vec_fast_precomp(poly, vals, tree, w, length, mod);
        prof_stop();

        _nmod_poly_tree_free(tree, length);
        flint_free(w);
    }

    _nmod_vec_clear(pts);
    _nmod_vec_clear(vals);
    _nmod_vec_clear(poly);

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
    slong length = info->length;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(pts, state, length, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            nn_ptr * tree = _nmod_poly_tree_alloc(length);
            _nmod_poly_tree_build(tree, pts, length, mod);

            nn_ptr w = _nmod_vec_init(length);
            _nmod_poly_interpolation_weights(w, tree, length, mod);

            _nmod_poly_tree_free(tree, length);
            flint_free(w);
        }
        prof_stop();
    }

    _nmod_vec_clear(pts);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, with precomp (by default, on `length` points) */
void sample_nmod_vec_geom_fast(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;

    FLINT_TEST_INIT(state);

    nn_ptr vals = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);
        nmod_poly_t poly;
        nmod_poly_init(poly, n);

        _nmod_vec_rand(vals, state, length, mod);
        ulong r = nmod_find_root(2*length, mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_poly_interpolate_geometric_nmod_vec_fast(poly, r, vals, length);
        prof_stop();
        nmod_poly_clear(poly);
    }

    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, not counting precomputations (which is on npoints_precomp points) */
void sample_nmod_vec_geom_fast_precomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong npoints_precomp = info->npoints_precomp;

    FLINT_TEST_INIT(state);

    nn_ptr vals = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        nn_ptr poly = _nmod_vec_init(length);
        _nmod_vec_rand(vals, state, length, mod);
        ulong r = nmod_find_root(2*npoints_precomp, mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        nmod_geometric_progression_t G;
        nmod_geometric_progression_init(G, r, npoints_precomp, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, vals, G, length, mod);
        prof_stop();

        nmod_geometric_progression_clear(G);
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, counting only precomputations (on npoints_precomp points) */
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
            _nmod_geometric_progression_init_function(G, r, npoints_precomp, mod, UWORD(2));
        prof_stop();

        nmod_geometric_progression_clear(G);
    }

    FLINT_TEST_CLEAR(state);
}

/* multiplying two polynomials */
void sample_nmod_poly_mul(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;

    FLINT_TEST_INIT(state);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        nmod_poly_t poly1;
        nmod_poly_t poly2;
        nmod_poly_t poly3;
        nmod_poly_init(poly1, n);
        nmod_poly_init(poly2, n);
        nmod_poly_init(poly3, n);
        nmod_poly_randtest_monic(poly1, state, length);
        nmod_poly_randtest_monic(poly2, state, length);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_poly_mul(poly3, poly1, poly2);
            /* nmod_poly_mulmid(poly3, poly1, poly2, length, 2*length); */
        prof_stop();

        nmod_poly_clear(poly1);
        nmod_poly_clear(poly2);
        nmod_poly_clear(poly3);
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
                     "              positive integer; precomputation will be on pre_fact * npoints\n"
                     "       [short] is optional (default 0)\n"
                     "               0 -> short bench, 1 -> full bench\n",
                     argv[0], argv[0]);
        return 0;
    }

    const slong func_bench = (argc >= 2) ? atoi(argv[1]) : 0;

    const slong pre_fact = (argc >= 3) ? atoi(argv[2]) : 1;

    if (pre_fact != 1)
    {
        flint_printf("[pre_fact provided] -> will build geometric progression with `%ld * points` points\n", pre_fact);
        flint_printf("[pre_fact provided] -> for general points (except precomp only), tree has `points` points\n");
    }

    /* number of lens differs: long test / short test */
    const slong nb_lens = (argc >= 4 && atoi(argv[3]) != 0) ? 25 : 21;
    slong lengths[] = {1, 2, 3, 4, 6,
                       8, 10, 12, 16, 20,
                       30, 45, 70, 100, 200,
                       400, 800, 1600, 3200, 6400,
                       12800, 25600, 51200, 102400, 204800};
    /* specific lengths for finding out fft_small threshold */
    /* slong lengths[] = { */
    /*     30, */
    /*     32, */
    /*     34, */
    /*     36, */
    /*     38, */
    /*     40, */
    /*     42, */
    /*     44, */
    /*     46, */
    /*     48, */
    /*     50, */
    /*     52, */
    /*     54, */
    /*     56, */
    /*     58, */
    /*     60, */
    /*     62, */
    /*     64, */
    /*     66, */
    /*     68, */
    /*     70, */
    /*     72, */
    /*     74, */
    /*     76, */
    /*     78, */
    /* }; */

    double tmp;
    double time_newton;
    double time_barycentric;
    double time_fast;
    double time_fast_precomp;
    double time_fast_onlyprecomp;
    double time_geom_fast;
    double time_geom_fast_precomp;
    double time_geom_fast_onlyprecomp;
    double time_poly_mul;

    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: measurements in ms\n");

    for (i = 63; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("==== nbits = %ld====\n", i);

        if (func_bench == 0)  /* bench all */
        {
            flint_printf("==== nb precomp points : %d * (poly length) ====\n", pre_fact);
            flint_printf("len\tpoints |            GENERAL POINTS              | GEOMETRIC PROGRESSION  | POLY_MUL\n");
            flint_printf("len\tpoints |newton\tbaryc\tfast\tw/ tree\ttree\t| fast\tw/ prec\tprecomp\n");
        }
        else if (func_bench == 1)  /* general only */
        {
            flint_printf("==== nb precomp points : %d * (poly length) ====\n", pre_fact);
            flint_printf("len\tpoints |            GENERAL POINTS              |\n");
            flint_printf("len\tpoints |newton\tbaryc\tfast\tw/ tree\ttree\t|\n");
        }
        else if (func_bench == 2)  /* geometric only */
        {
            flint_printf("==== nb precomp points : %d * (poly length) ====\n", pre_fact);
            flint_printf("len\tpoints | GEOMETRIC PROGRESSION  | POLY_MUL\n");
            flint_printf("len\tpoints |fast\tw/ prec\tprecomp\n");
        }

        for (int len = 0; len < nb_lens; ++len)
        {
            /* cycles / limb (if CLOCK_SCALE_FACTOR set correctly) */
            /* const double fac = npoints[len] * (double)FLINT_CLOCK_SCALE_FACTOR; */
            /* time in s */
            /* const double fac = 1000000. * __NB_ITER; */
            /* time in ms */
            const double fac = 1. * __NB_ITER;
            info.npoints_precomp = pre_fact * lengths[len];
            info.length = lengths[len];

            if (func_bench == 0 || func_bench == 1)
            {
                if (info.length <= 512)
                    prof_repeat(&time_newton, &tmp, sample_nmod_vec_newton, (void *) &info);
                else
                    time_newton = 0.0;

                if (info.length <= 1500)
                    prof_repeat(&time_barycentric, &tmp, sample_nmod_vec_barycentric, (void *) &info);
                else
                    time_barycentric = 0.0;

                prof_repeat(&time_fast, &tmp, sample_nmod_vec_fast, (void *) &info);

                prof_repeat(&time_fast_precomp, &tmp, sample_nmod_vec_fast_precomp, (void *) &info);

                prof_repeat(&time_fast_onlyprecomp, &tmp, sample_nmod_vec_fast_onlyprecomp, (void *) &info);
            }

            if (func_bench == 0 || func_bench == 2)
            {
                prof_repeat(&time_geom_fast, &tmp, sample_nmod_vec_geom_fast, (void *) &info);

                prof_repeat(&time_geom_fast_precomp, &tmp, sample_nmod_vec_geom_fast_precomp, (void *) &info);

                prof_repeat(&time_geom_fast_onlyprecomp, &tmp, sample_nmod_vec_geom_fast_onlyprecomp, (void *) &info);

                prof_repeat(&time_poly_mul, &tmp, sample_nmod_poly_mul, (void *) &info);
            }

            if (func_bench == 0)
            {
                flint_printf("%ld\t%7ld|%.1e %.1e %.1e %.1e %.1e |%.1e %.1e %.1e |%.1e\n",
                             info.length, info.npoints_precomp,
                             time_newton/fac, time_barycentric/fac, time_fast/fac, time_fast_precomp/fac, time_fast_onlyprecomp/fac,
                             time_geom_fast/fac, time_geom_fast_precomp/fac, time_geom_fast_onlyprecomp/fac,
                             time_poly_mul/fac);
            }
            else if (func_bench == 1)
            {
                flint_printf("%ld\t%7ld|%.1e %.1e %.1e %.1e %.1e\n",
                             info.length, info.npoints_precomp,
                             time_newton/fac, time_barycentric/fac, time_fast/fac, time_fast_precomp/fac, time_fast_onlyprecomp/fac);
            }
            else if (func_bench == 2)
            {
                flint_printf("%ld\t%7ld|%.1e %.1e %.1e |%.1e\n",
                             info.length, info.npoints_precomp,
                             time_geom_fast/fac, time_geom_fast_precomp/fac, time_geom_fast_onlyprecomp/fac,
                             time_poly_mul/fac);
            }
        }

        flint_printf("\n");
    }

    return 0;
}

#undef __NB_ITER
