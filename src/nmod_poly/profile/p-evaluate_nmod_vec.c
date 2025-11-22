/*
   Copyright (C) 2025 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* TODO profile when length << npoints */

#include "nmod.h"
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_vec.h"

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

        nn_ptr * tree = _nmod_poly_tree_alloc(npoints);
        _nmod_poly_tree_build(tree, pts, npoints, mod);

        prof_start();
        _nmod_poly_evaluate_nmod_vec_fast_precomp(vals, poly, length, tree, npoints, mod);
        prof_stop();

        _nmod_poly_tree_free(tree, npoints);
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
    slong npoints = info->npoints;

    FLINT_TEST_INIT(state);

    nn_ptr pts = _nmod_vec_init(npoints);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);
        _nmod_vec_rand(pts, state, npoints, mod);

        nn_ptr * tree = _nmod_poly_tree_alloc(npoints);

        prof_start();
        _nmod_poly_tree_build(tree, pts, npoints, mod);
        prof_stop();

        _nmod_poly_tree_free(tree, npoints);
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
        _nmod_poly_evaluate_geometric_nmod_vec_iter(vals, poly, length, r, npoints, mod);
        prof_stop();
        _nmod_vec_clear(poly);
    }

    _nmod_vec_clear(vals);

    FLINT_TEST_CLEAR(state);
}

/* geometric: fast, tree */
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

        nmod_geometric_progression_t G;
        nmod_geometric_progression_init(G, r, FLINT_MAX(npoints, length), mod);

        prof_start();
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
    slong length = info->length;
    slong npoints = info->npoints;

    FLINT_TEST_INIT(state);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        ulong r = nmod_find_root(2*npoints, mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        nmod_geometric_progression_t G;

        prof_start();
        nmod_geometric_progression_init(G, r, FLINT_MAX(npoints, length), mod);
        prof_stop();

        nmod_geometric_progression_clear(G);
    }

    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    const slong nb_lens = 25;
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

    flint_printf("unit: all measurements in c/l\n");

    for (i = 62; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("==== nbits = %ld====\n", i);
        flint_printf("len\tpoints |iter\tfast\tf-prcp\tf-tree |geom-it\tgeom-f\tgeom-pr\tgeom-precomp\n");
        for (int len = 0; len < nb_lens; ++len)
        {
            const double fac = lengths[len] * (double)FLINT_CLOCK_SCALE_FACTOR;
            info.length = lengths[len];
            info.npoints = lengths[len];

            if (info.npoints <= 1024)
            {
                prof_repeat(&time_iter, &tmp, sample_nmod_vec_iter, (void *) &info);
                time_iter /= fac;
            }
            else
                time_iter = 0.0;

            prof_repeat(&time_fast, &tmp, sample_nmod_vec_fast, (void *) &info);
            time_fast /= fac;

            prof_repeat(&time_fast_precomp, &tmp, sample_nmod_vec_fast_precomp, (void *) &info);
            time_fast_precomp /= fac;

            prof_repeat(&time_fast_onlyprecomp, &tmp, sample_nmod_vec_fast_onlyprecomp, (void *) &info);
            time_fast_onlyprecomp /= fac;

            if (info.npoints <= 1024)
            {
                prof_repeat(&time_geom_iter, &tmp, sample_nmod_vec_geom_iter, (void *) &info);
                time_geom_iter /= fac;
            }
            else
                time_geom_iter = 0.0;

            prof_repeat(&time_geom_fast, &tmp, sample_nmod_vec_geom_fast, (void *) &info);
            time_geom_fast /= fac;

            prof_repeat(&time_geom_fast_precomp, &tmp, sample_nmod_vec_geom_fast_precomp, (void *) &info);
            time_geom_fast_precomp /= fac;

            prof_repeat(&time_geom_fast_onlyprecomp, &tmp, sample_nmod_vec_geom_fast_onlyprecomp, (void *) &info);
            time_geom_fast_onlyprecomp /= fac;

            flint_printf("%ld\t%7ld|%.1e\t%.1e\t%.1e\t%.1e|%.1e\t%.1e\t%.1e\t%.1e\n",
                    info.length,
                    info.npoints,
                    time_iter,
                    time_fast,
                    time_fast_precomp,
                    time_fast_onlyprecomp,
                    time_geom_iter,
                    time_geom_fast,
                    time_geom_fast_precomp,
                    time_geom_fast_onlyprecomp);
        }

        flint_printf("\n");
    }

    return 0;
}
