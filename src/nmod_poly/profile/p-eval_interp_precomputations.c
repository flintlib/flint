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
static ulong nmod_find_root(slong n, nmod_t mod)
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

/* precomputation for general points (subproduct tree) */
void sample_general_precomp(void * arg, ulong count)
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

/* precomputations for geometric progression (eval+interp+extrap) */
void sample_geometric_precomp(void * arg, ulong count)
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
        {
            nmod_geometric_progression_init(G, r, npoints_precomp, mod);
            nmod_geometric_progression_clear(G);
        }
        prof_stop();

    }

    FLINT_TEST_CLEAR(state);
}

/* precomputations for geometric progression (eval only) */
void sample_geometric_precomp_eval(void * arg, ulong count)
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
        G->len = npoints_precomp;
        G->mod = mod;

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_geometric_progression_init_function(G, r, npoints_precomp, mod, UWORD(1));
            nmod_geometric_progression_clear(G);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/* precomputations for geometric progression (interp only) */
void sample_geometric_precomp_interp(void * arg, ulong count)
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
        G->len = npoints_precomp;
        G->mod = mod;

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_geometric_progression_init_function(G, r, npoints_precomp, mod, UWORD(2));
            nmod_geometric_progression_clear(G);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/* precomputations for geometric progression (extrap only) */
void sample_geometric_precomp_extrap(void * arg, ulong count)
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
        G->len = npoints_precomp;
        G->mod = mod;

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_geometric_progression_init_function(G, r, npoints_precomp, mod, UWORD(4));
            nmod_geometric_progression_clear(G);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/* precomputations for geometric progression (extra+inter only) */
void sample_geometric_precomp_extra_inter(void * arg, ulong count)
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
        G->len = npoints_precomp;
        G->mod = mod;

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_geometric_progression_init_function(G, r, npoints_precomp, mod, UWORD(6));
            nmod_geometric_progression_clear(G);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

int main(int argc, char * argv[])
{
    if (argc > 1 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        flint_printf("Usage: %s -h for this help, or\n"
                     "       %s [func] [short]\n"
                     " Optional arguments (if one is provided, previous ones must as well)\n"
                     "       [func] is optional (default 0)\n"
                     "              0 -> all\n"
                     "              1 -> general points only\n"
                     "              2 -> geometric points only\n"
                     "       [short] is optional (default 0)\n"
                     "               0 -> short bench, 1 -> full bench\n",
                     argv[0], argv[0]);
        return 0;
    }

    const slong func_bench = (argc >= 2) ? atoi(argv[1]) : 0;

    /* number of lens differs: long test / short test */
    const slong nb_lens = (argc >= 3 && atoi(argv[2]) != 0) ? 25 : 21;
    slong lengths[] = {1, 2, 3, 4, 6,
                       8, 10, 12, 16, 20,
                       30, 45, 70, 100, 200,
                       400, 800, 1600, 3200, 6400,
                       12800, 25600, 51200, 102400, 204800};

    double tmp;
    double time_general;
    double time_geometric_all;
    double time_geometric_eval;
    double time_geometric_interp;
    double time_geometric_extrap;
    double time_geometric_extra_inter;

    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: measurements in ms\n");

    for (i = 63; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("==== nbits = %ld====\n", i);

        if (func_bench == 0)  /* bench all */
            flint_printf("len\tpoints | general |   geom  |  g-eval | g-interp| g-extrap| g-int-ext \n");
        else if (func_bench == 1)  /* general only */
            flint_printf("len\tpoints | general\n");
        else if (func_bench == 2)  /* geometric only */
            flint_printf("len\tpoints |   geom  |  g-eval | g-interp| g-extrap| g-int-ext \n");

        for (int len = 0; len < nb_lens; ++len)
        {
            /* time in microsec */
            const double fac = 1. * __NB_ITER;
            info.npoints_precomp = lengths[len];
            info.length = lengths[len];

            if (func_bench == 0 || func_bench == 1)
                prof_repeat(&time_general, &tmp, sample_general_precomp, (void *) &info);

            if (func_bench == 0 || func_bench == 2)
            {
                prof_repeat(&time_geometric_all, &tmp, sample_geometric_precomp, (void *) &info);
                prof_repeat(&time_geometric_eval, &tmp, sample_geometric_precomp_eval, (void *) &info);
                prof_repeat(&time_geometric_interp, &tmp, sample_geometric_precomp_interp, (void *) &info);
                prof_repeat(&time_geometric_extrap, &tmp, sample_geometric_precomp_extrap, (void *) &info);
                prof_repeat(&time_geometric_extra_inter, &tmp, sample_geometric_precomp_extra_inter, (void *) &info);
            }

            if (func_bench == 0)
            {
                flint_printf("%ld\t%7ld| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e\n",
                             info.length, info.npoints_precomp,
                             time_general/fac,
                             time_geometric_all/fac,
                             time_geometric_eval/fac,
                             time_geometric_interp/fac,
                             time_geometric_extrap/fac,
                             time_geometric_extra_inter/fac);
            }
            else if (func_bench == 1)
            {
                flint_printf("%ld\t%7ld| %.1e\n",
                             info.length, info.npoints_precomp,
                             time_general/fac);
            }
            else if (func_bench == 2)
            {
                flint_printf("%ld\t%7ld| %.1e | %.1e | %.1e | %.1e | %.1e\n",
                             info.length, info.npoints_precomp,
                             time_geometric_all/fac,
                             time_geometric_eval/fac,
                             time_geometric_interp/fac,
                             time_geometric_extrap/fac,
                             time_geometric_extra_inter/fac);
            }
        }

        flint_printf("\n");
    }

    return 0;
}

#undef __NB_ITER
