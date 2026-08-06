/*
   Copyright (C) 2026 Vincent Neiger, Kevin Tran

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
    slong ilen;
    slong olen;
} info_t;

/* general points */
void sample_general(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ipts = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr opts = FLINT_ARRAY_ALLOC(olen, ulong);
    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);

        /* need pairwise distinct ipts for interpolation, let's take ipts[k] = k */
        for (slong k = 0; k < ilen; k++)
            ipts[k] = k;
        _nmod_vec_rand(opts, state, olen, mod);
        _nmod_vec_rand(ival, state, ilen, mod);

        prof_start();
        nn_ptr poly = FLINT_ARRAY_ALLOC(ilen, ulong);
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_poly_interpolate_nmod_vec(poly, ipts, ival, ilen, mod);
            _nmod_poly_evaluate_nmod_vec(oval, poly, ilen, opts, olen, mod);
        }
        flint_free(poly);
        prof_stop();
    }

    flint_free(ipts);
    flint_free(opts);
    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* general points, not counting tree precomputations */
void sample_general_noprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ipts = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr opts = FLINT_ARRAY_ALLOC(olen, ulong);
    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);

        /* need pairwise distinct ipts for interpolation, let's take ipts[k] = k */
        for (slong k = 0; k < ilen; k++)
            ipts[k] = k;
        _nmod_vec_rand(opts, state, olen, mod);
        _nmod_vec_rand(ival, state, ilen, mod);

        nn_ptr * itree = _nmod_poly_tree_alloc(ilen);
        _nmod_poly_tree_build(itree, ipts, ilen, mod);
        nn_ptr iw = _nmod_vec_init(ilen);
        _nmod_poly_interpolation_weights(iw, itree, ilen, mod);

        nn_ptr * otree = _nmod_poly_tree_alloc(olen);
        _nmod_poly_tree_build(otree, opts, olen, mod);
        nn_ptr ow = _nmod_vec_init(olen);
        _nmod_poly_interpolation_weights(ow, otree, olen, mod);

        prof_start();
        nn_ptr poly = FLINT_ARRAY_ALLOC(ilen, ulong);
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_poly_interpolate_nmod_vec_fast_precomp(poly, ival, itree, iw, ilen, mod);
            _nmod_poly_evaluate_nmod_vec_fast_precomp(oval, poly, ilen, otree, olen, mod);
        }
        flint_free(poly);
        prof_stop();

        _nmod_poly_tree_free(itree, ilen);
        flint_free(iw);
        _nmod_poly_tree_free(otree, olen);
        flint_free(ow);
    }

    flint_free(ipts);
    flint_free(opts);
    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* tree, counting only precomputations */
void sample_general_onlyprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ipts = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr opts = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        if (n == UWORD(0)) n++;
        nmod_init(&mod, n);

        /* need pairwise distinct ipts for interpolation, let's take ipts[k] = k */
        for (slong k = 0; k < ilen; k++)
            ipts[k] = k;
        _nmod_vec_rand(opts, state, olen, mod);

        nn_ptr * itree = _nmod_poly_tree_alloc(ilen);
        nn_ptr iw = _nmod_vec_init(ilen);
        nn_ptr * otree = _nmod_poly_tree_alloc(olen);
        nn_ptr ow = _nmod_vec_init(olen);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_poly_tree_build(itree, ipts, ilen, mod);
            _nmod_poly_interpolation_weights(iw, itree, ilen, mod);

            _nmod_poly_tree_build(otree, opts, olen, mod);
            _nmod_poly_interpolation_weights(ow, otree, olen, mod);
        }
        prof_stop();
        _nmod_poly_tree_free(itree, ilen);
        flint_free(iw);
        _nmod_poly_tree_free(otree, olen);
        flint_free(ow);
    }

    flint_free(ipts);
    flint_free(opts);

    FLINT_TEST_CLEAR(state);
}

/* geometric, including precomp */
void sample_geometric(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);

        _nmod_vec_rand(ival, state, ilen, mod);
        ulong r = nmod_find_root(2*(ilen + olen), mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_poly_extrapolate_geometric(oval, olen, ival, ilen, ilen, r, mod);
        prof_stop();
    }

    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* geometric, not counting precomputations */
void sample_geometric_noprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        nmod_geometric_progression_t G;
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);

        _nmod_vec_rand(ival, state, ilen, mod);
        ulong r = nmod_find_root(2*(ilen + olen), mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");
        _nmod_geometric_progression_init_function(G, r, ilen+olen, mod, UWORD(4));

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_poly_extrapolate_geometric_precomp(oval, olen, ival, ilen, ilen, G);
        prof_stop();

        nmod_geometric_progression_clear(G);
    }

    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* geometric (backward), including precomp */
void sample_geometric_backward(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);

        _nmod_vec_rand(ival, state, ilen, mod);
        ulong r = nmod_find_root(2*(ilen + olen), mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_poly_extrapolate_geometric(oval, olen, ival, ilen, -olen, r, mod);
        prof_stop();
    }

    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* geometric, not counting precomputations */
void sample_geometric_backward_noprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        nmod_geometric_progression_t G;
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);

        _nmod_vec_rand(ival, state, ilen, mod);
        ulong r = nmod_find_root(2*(ilen + olen), mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");
        _nmod_geometric_progression_init_function(G, r, ilen+olen, mod, UWORD(4));

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            nmod_poly_extrapolate_geometric_precomp(oval, olen, ival, ilen, -olen, G);
        prof_stop();

        nmod_geometric_progression_clear(G);
    }

    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* geometric, not counting precomputations, via interpolation+evaluation */
void sample_geometric_noprecomp_intev(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    nn_ptr ival = FLINT_ARRAY_ALLOC(ilen, ulong);
    nn_ptr oval = FLINT_ARRAY_ALLOC(olen, ulong);

    for (i = 0; i < count; i++)
    {
        nmod_geometric_progression_t G;
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);

        _nmod_vec_rand(ival, state, ilen, mod);
        ulong r = nmod_find_root(2*(ilen + olen), mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");
        _nmod_geometric_progression_init_function(G, r, ilen+olen, mod, UWORD(3));

        prof_start();
        nn_ptr poly = FLINT_ARRAY_ALLOC(ilen, ulong);
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, ival, G, ilen, mod);
            /* scale so that next evaluations are at q**j for j=ilen...ilen+olen-1 */
            ulong q_pow_ilen = nmod_pow_ui(r, 2*ilen, mod);
            ulong qk = q_pow_ilen;
            for (slong k = 1; k < ilen; k++)
            {
                poly[k] = nmod_mul(poly[k], qk, mod);
                qk = nmod_mul(qk, q_pow_ilen, mod);
            }
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(oval, poly, ilen, G, olen, mod);
        }
        flint_free(poly);
        prof_stop();

        nmod_geometric_progression_clear(G);
    }

    flint_free(ival);
    flint_free(oval);

    FLINT_TEST_CLEAR(state);
}

/* geometric, counting only precomputations (on npoints_precomp points) */
void sample_geometric_onlyprecomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    for (i = 0; i < count; i++)
    {
        nmod_geometric_progression_t G;
        n = n_randprime(state, bits, 0);
        nmod_init(&mod, n);

        ulong r = nmod_find_root(2*(ilen + olen), mod);
        if (r == 0)
            flint_printf("\n...could not find element of suitable order for geometric progression...\n");

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
        {
            _nmod_geometric_progression_init_function(G, r, ilen+olen, mod, UWORD(4));
            nmod_geometric_progression_clear(G);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/* comparison with relevant middle product */
void sample_nmod_poly_mulmid(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    ulong i;

    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong ilen = info->ilen;
    slong olen = info->olen;

    FLINT_TEST_INIT(state);

    for (i = 0; i < count; i++)
    {
        n = n_randprime(state, bits, 1);
        nmod_init(&mod, n);
        nn_ptr poly1 = FLINT_ARRAY_ALLOC(ilen+olen-1, ulong);
        nn_ptr poly2 = FLINT_ARRAY_ALLOC(ilen, ulong);
        nn_ptr poly3 = FLINT_ARRAY_ALLOC(olen, ulong);

        _nmod_vec_rand(poly1, state, ilen+olen-1, mod);
        _nmod_vec_rand(poly2, state, ilen, mod);

        prof_start();
        for (ulong ii = 0; ii < __NB_ITER; ii++)
            _nmod_poly_mulmid(poly3, poly1, ilen+olen-1, poly2, ilen, ilen-1, ilen+olen-1, mod);
        prof_stop();

        flint_free(poly1);
        flint_free(poly2);
        flint_free(poly3);
    }

    FLINT_TEST_CLEAR(state);
}

int main(int argc, char * argv[])
{
    if (argc > 1 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        flint_printf("Usage: %s -h for this help, or\n"
                     "       %s [func] [growth] [short]\n"
                     " Optional arguments (if one is provided, previous ones must as well)\n"
                     "       [func] is optional (default 0)\n"
                     "              0 -> all\n"
                     "              1 -> general points only\n"
                     "              2 -> geometric points only\n"
                     "       [growth] is optional (default 1.0)\n"
                     "              positive float: extrapolate from `ilen` points to `growth * ilen` points\n"
                     "       [short] is optional (default 0)\n"
                     "               0 -> short bench, 1 -> full bench\n",
                     argv[0], argv[0]);
        return 0;
    }

    const slong func_bench = (argc >= 2) ? atoi(argv[1]) : 0;

    const double growth = (argc >= 3) ? atof(argv[2]) : 1.;

    if (growth != 1.)
        flint_printf("[growth provided] -> will extrapolate from `ilen` points to `%lf * ilen` points\n", growth);

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
    double time_general;
    double time_general_noprecomp;
    double time_general_onlyprecomp;
    double time_geometric;
    double time_geometric_noprecomp;
    double time_geometric_backward;
    double time_geometric_backward_noprecomp;
    double time_geometric_intev;
    double time_geometric_onlyprecomp;
    double time_poly_mulmid;

    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: measurements in ms\n");

    for (i = 63; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("==== nbits = %ld====\n", i);

        if (func_bench == 0)  /* bench all */
            flint_printf("ilen\t  olen |general\tw/ tree\t  tree | geom\tprec\tgeom-b\tprec-b\tint+ev |precomp\t|mulmid\n");
        else if (func_bench == 1)  /* general only */
            flint_printf("ilen\t  olen |general\tw/ tree\t  tree |\n");
        else if (func_bench == 2)  /* geometric only */
            flint_printf("ilen\t  olen |geom\tprec\tgeom-b\tprec-b\tint+ev |precomp\t|mulmid\n");

        for (int len = 0; len < nb_lens; ++len)
        {
            /* cycles / limb (if CLOCK_SCALE_FACTOR set correctly) */
            /* const double fac = npoints[len] * (double)FLINT_CLOCK_SCALE_FACTOR; */
            /* time in s */
            /* const double fac = 1000000. * __NB_ITER; */
            /* time in ms */
            const double fac = 1. * __NB_ITER;
            info.ilen = lengths[len];
            info.olen = growth * lengths[len];

            if (func_bench == 0 || func_bench == 1)
            {
                prof_repeat(&time_general, &tmp, sample_general, (void *) &info);

                prof_repeat(&time_general_noprecomp, &tmp, sample_general_noprecomp, (void *) &info);

                prof_repeat(&time_general_onlyprecomp, &tmp, sample_general_onlyprecomp, (void *) &info);
            }

            if (func_bench == 0 || func_bench == 2)
            {
                prof_repeat(&time_geometric, &tmp, sample_geometric, (void *) &info);

                prof_repeat(&time_geometric_noprecomp, &tmp, sample_geometric_noprecomp, (void *) &info);

                prof_repeat(&time_geometric_backward, &tmp, sample_geometric_backward, (void *) &info);

                prof_repeat(&time_geometric_backward_noprecomp, &tmp, sample_geometric_backward_noprecomp, (void *) &info);

                prof_repeat(&time_geometric_intev, &tmp, sample_geometric_noprecomp_intev, (void *) &info);

                prof_repeat(&time_geometric_onlyprecomp, &tmp, sample_geometric_onlyprecomp, (void *) &info);

                prof_repeat(&time_poly_mulmid, &tmp, sample_nmod_poly_mulmid, (void *) &info);
            }

            if (func_bench == 0)
            {
                flint_printf("%ld\t%7ld|%.1e %.1e %.1e|%.1e %.1e %.1e %.1e %.1e|%.1e |%.1e\n",
                             info.ilen, info.olen,
                             time_general/fac, time_general_noprecomp/fac, time_general_onlyprecomp/fac,
                             time_geometric/fac, time_geometric_noprecomp/fac,
                             time_geometric_backward/fac, time_geometric_backward_noprecomp/fac,
                             time_geometric_intev/fac, time_geometric_onlyprecomp/fac,
                             time_poly_mulmid/fac);
            }
            else if (func_bench == 1)
            {
                flint_printf("%ld\t%7ld|%.1e %.1e %.1e\n",
                             info.ilen, info.olen,
                             time_general/fac, time_general_noprecomp/fac, time_general_onlyprecomp/fac);
            }
            else if (func_bench == 2)
            {
                flint_printf("%ld\t%7ld|%.1e %.1e %.1e %.1e %.1e|%.1e |%.1e\n",
                             info.ilen, info.olen,
                             time_geometric/fac, time_geometric_noprecomp/fac,
                             time_geometric_backward/fac, time_geometric_backward_noprecomp/fac,
                             time_geometric_intev/fac, time_geometric_onlyprecomp/fac,
                             time_poly_mulmid/fac);
            }
        }

        flint_printf("\n");
    }

    return 0;
}

#undef __NB_ITER
