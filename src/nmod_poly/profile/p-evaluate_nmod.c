/*
   Copyright (C) 2024 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_vec.h"

#define NB_LENS 20
#define NB_REPS 10

/* only set if wanting the full bench with static functions */
/* --> this implies temporarily making non-static in the relevant files */
#define __DO_BENCH_STATIC 0
#if __DO_BENCH_STATIC  /* declare prototypes */
ulong _nmod_poly_evaluate_one(nn_srcptr poly, slong len, ulong modn);
ulong _nmod_poly_evaluate_one1(nn_srcptr poly, slong len, ulong modn);
ulong _nmod_poly_evaluate_one2(nn_srcptr poly, slong len, ulong modn);
ulong _nmod_poly_evaluate_mone(nn_srcptr poly, slong len, ulong modn);
ulong _nmod_poly_evaluate_mone1(nn_srcptr poly, slong len, ulong modn);
ulong _nmod_poly_evaluate_mone2(nn_srcptr poly, slong len, ulong modn);
#endif

typedef struct
{
    flint_bitcnt_t bits;
    slong length;
    slong ctype;  /* -1: c == -1; +1: c == +1; else: c random */
} info_t;

void sample_interface(void * arg, ulong count)
{
    ulong n;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong ctype = info->ctype;
    slong i, j;

    nmod_poly_t poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    for (i = 0; i < (slong)count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_poly_init2(poly, n, length);
        _nmod_poly_set_length(poly, length);
        for (j = 0; j < length; j++)
            poly->coeffs[j] = n_randint(state, n);

        if (ctype == -1)
            pt = n-1;
        else if (ctype == 1)
            pt = UWORD(1);
        else
            pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < NB_REPS; j++)
            nmod_poly_evaluate_nmod(poly, pt);
        prof_stop();

        nmod_poly_clear(poly);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_generic(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong ctype = info->ctype;
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < (slong)count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        if (ctype == -1)
            pt = n-1;
        else if (ctype == 1)
            pt = UWORD(1);
        else
            pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < NB_REPS; j++)
            _nmod_poly_evaluate_nmod(poly, length, pt, mod);
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

void sample_precomp(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong ctype = info->ctype;
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < (slong)count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        if (ctype == -1)
            pt = n-1;
        else if (ctype == 1)
            pt = UWORD(1);
        else
            pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < NB_REPS; j++)
        {
            const ulong pt_precomp = n_mulmod_precomp_shoup(pt, n);
            _nmod_poly_evaluate_nmod_precomp(poly, length, pt, pt_precomp, n);
        }
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

void sample_precomp_lazy(void * arg, ulong count)
{
    ulong n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong ctype = info->ctype;
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < (slong)count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        if (ctype == -1)
            pt = n-1;
        else if (ctype == 1)
            pt = UWORD(1);
        else
            pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < NB_REPS; j++)
        {
            const ulong pt_precomp = n_mulmod_precomp_shoup(pt, n);
            ulong val = _nmod_poly_evaluate_nmod_precomp_lazy(poly, length, pt, pt_precomp, n);
            if (val >= 2*n)
                val -= 2*n;
            else if (val >= n)
                val -= n;
        }
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

#if __DO_BENCH_STATIC

#define SAMPLE_ONE_MONE(variant)                \
void sample_##variant(void * arg, ulong count)  \
{                                               \
    ulong n;                                    \
    info_t * info = (info_t *) arg;             \
    flint_bitcnt_t bits = info->bits;           \
    slong length = info->length;                \
    slong i, j;                                 \
                                                \
    ulong * poly;                               \
                                                \
    FLINT_TEST_INIT(state);                     \
                                                \
    poly = _nmod_vec_init(length);              \
                                                \
    for (i = 0; i < (slong)count; i++)          \
    {                                           \
        n = n_randbits(state, bits);            \
        if (n == UWORD(0)) n++;                 \
                                                \
        for (j = 0; j < length; j++)            \
            poly[j] = n_randint(state, n);      \
                                                \
        prof_start();                           \
        for (j = 0; j < NB_REPS; j++)           \
            _nmod_poly_evaluate_##variant(      \
                        poly, length, n);       \
        prof_stop();                            \
    }                                           \
                                                \
    _nmod_vec_clear(poly);                      \
                                                \
    FLINT_TEST_CLEAR(state);                    \
}

SAMPLE_ONE_MONE(one)
SAMPLE_ONE_MONE(one1)
SAMPLE_ONE_MONE(one2)
SAMPLE_ONE_MONE(mone)
SAMPLE_ONE_MONE(mone1)
SAMPLE_ONE_MONE(mone2)
#endif

int main(void)
{
    slong lengths[NB_LENS] = 
        { 1, 2, 3, 4, 6, 8,
          10, 12, 16, 20, 32, 45,
          64, 128, 256, 1024,
          8192, 65536, 200000, 1000000};

    double min, max;
    double mins_interface[NB_LENS];
    double mins_interface_one[NB_LENS];
    double mins_interface_mone[NB_LENS];
    double mins_generic[NB_LENS];
    double mins_precomp[NB_LENS];
    double mins_precomp_lazy[NB_LENS];
#if __DO_BENCH_STATIC
    double mins_one[NB_LENS];
    double mins_one1[NB_LENS];
    double mins_one2[NB_LENS];
    double mins_mone[NB_LENS];
    double mins_mone1[NB_LENS];
    double mins_mone2[NB_LENS];
#endif
    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: all measurements in c/l\n");
    flint_printf("profiled: interface | generic | precomp | precomp_lazy\n");

    for (i = 62; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("nbits = %ld\n", i);
        for (int len = 0; len < NB_LENS; ++len)
        {
            info.length = lengths[len];

            info.ctype = 0;
            prof_repeat(&min, &max, sample_interface, (void *) &info);
            mins_interface[len-1] = min;
            prof_repeat(&min, &max, sample_generic, (void *) &info);
            mins_generic[len-1] = min;
            prof_repeat(&min, &max, sample_precomp, (void *) &info);
            mins_precomp[len-1] = min;
            prof_repeat(&min, &max, sample_precomp_lazy, (void *) &info);
            mins_precomp_lazy[len-1] = min;

            info.ctype = 1;
            prof_repeat(&min, &max, sample_interface, (void *) &info);
            mins_interface_one[len-1] = min;

            info.ctype = -1;
            prof_repeat(&min, &max, sample_interface, (void *) &info);
            mins_interface_mone[len-1] = min;

#if __DO_BENCH_STATIC
            prof_repeat(&min, &max, sample_one, (void *) &info);
            mins_one[len-1] = min;
            prof_repeat(&min, &max, sample_one1, (void *) &info);
            mins_one1[len-1] = min;
            prof_repeat(&min, &max, sample_one2, (void *) &info);
            mins_one2[len-1] = min;
            prof_repeat(&min, &max, sample_mone, (void *) &info);
            mins_mone[len-1] = min;
            prof_repeat(&min, &max, sample_mone1, (void *) &info);
            mins_mone1[len-1] = min;
            prof_repeat(&min, &max, sample_mone2, (void *) &info);
            mins_mone2[len-1] = min;
#endif
        }

#if __DO_BENCH_STATIC
        if (i < FLINT_BITS-1)
        {
            for (int len = 0; len < NB_LENS; ++len)
            {
                double fac = NB_REPS * (double)lengths[len] * (double)FLINT_CLOCK_SCALE_FACTOR;
                flint_printf("   len %ld\t\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf",
                        lengths[len],
                        mins_interface[len-1]/fac,
                        mins_generic[len-1]/fac,
                        mins_precomp[len-1]/fac,
                        mins_precomp_lazy[len-1]/fac,
                        mins_one[len-1]/fac,
                        mins_one1[len-1]/fac,
                        mins_one2[len-1]/fac,
                        mins_mone[len-1]/fac,
                        mins_mone1[len-1]/fac,
                        mins_mone2[len-1]/fac);
                flint_printf("\n");
            }
        }
#else
        if (i < FLINT_BITS-1)
        {
            flint_printf("                interf  generic precomp lazy    one     mone\n");
            for (int len = 0; len < NB_LENS; ++len)
            {
                double fac = NB_REPS * (double)lengths[len] * (double)FLINT_CLOCK_SCALE_FACTOR;
                flint_printf("   len %ld\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf",
                        lengths[len],
                        mins_interface[len-1]/fac,
                        mins_generic[len-1]/fac,
                        mins_precomp[len-1]/fac,
                        mins_precomp_lazy[len-1]/fac,
                        mins_interface_one[len-1]/fac,
                        mins_interface_mone[len-1]/fac);
                flint_printf("\n");
            }
        }
#endif
        if (i == FLINT_BITS-1)
        {
            flint_printf("                interf  generic precomp one     mone\n");
            for (int len = 0; len < NB_LENS; ++len)
            {
                double fac = NB_REPS * (double)lengths[len] * (double)FLINT_CLOCK_SCALE_FACTOR;
                flint_printf("   len %ld\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf",
                        lengths[len],
                        mins_interface[len-1]/fac,
                        mins_generic[len-1]/fac,
                        mins_precomp[len-1]/fac,
                        mins_interface_one[len-1]/fac,
                        mins_interface_mone[len-1]/fac);
                flint_printf("\n");
            }
        }

        if (i == FLINT_BITS)
        {
            flint_printf("                interf  generic one     mone\n");
            for (int len = 0; len < NB_LENS; ++len)
            {
                double fac = NB_REPS * (double)lengths[len] * (double)FLINT_CLOCK_SCALE_FACTOR;
                flint_printf("   len %ld\t%.2lf\t%.2lf\t%.2lf\t%.2lf",
                        lengths[len],
                        mins_interface[len-1]/fac,
                        mins_generic[len-1]/fac,
                        mins_interface_one[len-1]/fac,
                        mins_interface_mone[len-1]/fac);
                flint_printf("\n");
            }
        }

        flint_printf("\n");
    }

    return 0;
}

#undef __DO_BENCH_STATIC
