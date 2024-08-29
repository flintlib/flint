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

typedef struct
{
    flint_bitcnt_t bits;
    slong length;
} info_t;

void sample_interface(void * arg, ulong count)
{
    ulong n;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    slong i, j;

    nmod_poly_t poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_poly_init2(poly, n, length);
        _nmod_poly_set_length(poly, length);
        for (j = 0; j < length; j++)
            poly->coeffs[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
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
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
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
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
        {
            const ulong pt_precomp = n_mulmod_precomp_shoup(pt, n);
            _nmod_poly_evaluate_nmod_precomp(poly, length, pt, pt_precomp, mod);
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
    slong i, j;

    ulong * poly;
    ulong pt;

    FLINT_TEST_INIT(state);

    poly = _nmod_vec_init(length);

    for (i = 0; i < count; i++)
    {
        n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        nmod_init(&mod, n);

        for (j = 0; j < length; j++)
            poly[j] = n_randint(state, n);

        pt = n_randint(state, n);

        prof_start();
        for (j = 0; j < 100; j++)
        {
            const ulong pt_precomp = n_mulmod_precomp_shoup(pt, n);
            _nmod_poly_evaluate_nmod_precomp_lazy(poly, length, pt, pt_precomp, mod);
        }
        prof_stop();
    }

    _nmod_vec_clear(poly);

    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    slong lengths[18] = {1, 2, 3, 4, 6, 8,
                         10, 12, 16, 20, 32, 45,
                         64, 128, 256, 1024, 8192, 65536};

    double min, max;
    double mins[18]; // note: max seems to be consistently identical or extremely close to min
    double mins_generic[18];
    double mins_precomp[18];
    double mins_precomp_lazy[18];
    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: all measurements in c/l\n");
    flint_printf("profiled: interface | generic | precomp | precomp_lazy\n");

    for (i = 62; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        printf("nbits = %ld\n", i);
        for (int len = 0; len < 18; ++len)
        {
            info.length = lengths[len];

            prof_repeat(&min, &max, sample_interface, (void *) &info);
            mins[len-1] = min;
            prof_repeat(&min, &max, sample_generic, (void *) &info);
            mins_generic[len-1] = min;
            prof_repeat(&min, &max, sample_precomp, (void *) &info);
            mins_precomp[len-1] = min;
            prof_repeat(&min, &max, sample_precomp_lazy, (void *) &info);
            mins_precomp_lazy[len-1] = min;
        }

        if (i < FLINT_BITS-1)
        {
            for (int len = 0; len < 18; ++len)
            {
                flint_printf("   len %ld\t%.2lf\t%.2lf\t%.2lf\t%.2lf",
                        lengths[len],
                        (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_precomp[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_precomp_lazy[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100));
                flint_printf("\n");
            }
        }
        else if (i < FLINT_BITS)
        {
            for (int len = 0; len < 18; ++len)
            {
                flint_printf("   len %ld\t%.2lf\t%.2lf\t%.2lf\t na",
                        lengths[len],
                        (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_precomp[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_precomp_lazy[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100));
                flint_printf("\n");
            }
        }
        else  // i == FLINT_BITS
        {
            for (int len = 0; len < 18; ++len)
            {
                flint_printf("   len %ld\t%.2lf\t%.2lf\t na \t na",
                        lengths[len],
                        (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100),
                        (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(lengths[len]*100));
                flint_printf("\n");
            }
        }

        flint_printf("\n");
    }

    return 0;
}
