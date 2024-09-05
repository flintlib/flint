/*
   Copyright (C) 2024 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
   */

#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

void _nmod_mat_scalar_addmul_ui_generic(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B,
                                       const ulong c);
void _nmod_mat_scalar_addmul_ui_precomp(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B,
                                       const ulong c, const ulong c_pr);

typedef struct
{
    flint_bitcnt_t bits;
    slong length;
} info_t;

void sample(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    nmod_mat_t mat1, mat2, res;

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        ulong c = n_randint(state, n);

        nmod_mat_init(mat1, length, length, n);
        nmod_mat_init(mat2, length, length, n);
        nmod_mat_randtest(mat1, state);
        nmod_mat_randtest(mat2, state);

        nmod_mat_init(res, length, length, n);

        prof_start();
        for (slong j = 0; j < 100; j++)
            nmod_mat_scalar_addmul_ui(res, mat1, mat2, c);
        prof_stop();

        nmod_mat_clear(mat1);
        nmod_mat_clear(mat2);
        nmod_mat_clear(res);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_generic(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    nmod_mat_t mat1, mat2, res;

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        ulong c = n_randint(state, n);

        nmod_mat_init(mat1, length, length, n);
        nmod_mat_init(mat2, length, length, n);
        nmod_mat_randtest(mat1, state);
        nmod_mat_randtest(mat2, state);

        nmod_mat_init(res, length, length, n);

        prof_start();
        for (slong j = 0; j < 100; j++)
            _nmod_mat_scalar_addmul_ui_generic(res, mat1, mat2, c);
        prof_stop();

        nmod_mat_clear(mat1);
        nmod_mat_clear(mat2);
        nmod_mat_clear(res);
    }

    FLINT_TEST_CLEAR(state);
}

void sample_precomp(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    slong length = info->length;
    nmod_mat_t mat1, mat2, res;

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong n = n_randbits(state, bits);
        if (n == UWORD(0)) n++;

        ulong c = n_randint(state, n);

        nmod_mat_init(mat1, length, length, n);
        nmod_mat_init(mat2, length, length, n);
        nmod_mat_randtest(mat1, state);
        nmod_mat_randtest(mat2, state);

        nmod_mat_init(res, length, length, n);

        prof_start();
        for (slong j = 0; j < 100; j++)
        {
            const ulong c_pr = n_mulmod_precomp_shoup(c, n);
            _nmod_mat_scalar_addmul_ui_precomp(res, mat1, mat2, c, c_pr);
        }
        prof_stop();

        nmod_mat_clear(mat1);
        nmod_mat_clear(mat2);
        nmod_mat_clear(res);
    }

    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;
    double mins[18]; // note: max seems to be consistently identical or extremely close to min
    double mins_precomp[18]; // note: max seems to be consistently identical or extremely close to min
    double mins_generic[18]; // note: max seems to be consistently identical or extremely close to min
    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: all measurements in c/l\n");
    flint_printf("profiled: general interface | precomp shoup | generic\n");
    flint_printf("bit/len\t");
    for (int len = 1; len <= 16; ++len)
        flint_printf("%d\t\t", len);
    flint_printf("32\t\t");
    flint_printf("256\n");

    for (i = 2; i <= FLINT_BITS; i++)
    {
        info.bits = i;

        for (int len = 1; len <= 16; ++len)
        {
            info.length = len;

            prof_repeat(&min, &max, sample, (void *) &info);
            mins[len-1] = min;
            prof_repeat(&min, &max, sample_precomp, (void *) &info);
            mins_precomp[len-1] = min;
            prof_repeat(&min, &max, sample_generic, (void *) &info);
            mins_generic[len-1] = min;
        }

        info.length = 32;
        prof_repeat(&min, &max, sample, (void *) &info);
        mins[16] = min;
        prof_repeat(&min, &max, sample_precomp, (void *) &info);
        mins_precomp[16] = min;
        prof_repeat(&min, &max, sample_generic, (void *) &info);
        mins_generic[16] = min;

        info.length = 256;
        prof_repeat(&min, &max, sample, (void *) &info);
        mins[17] = min;
        prof_repeat(&min, &max, sample_precomp, (void *) &info);
        mins_precomp[17] = min;
        prof_repeat(&min, &max, sample_generic, (void *) &info);
        mins_generic[17] = min;

        if (i < FLINT_BITS)
        {
            flint_printf("%wd", i);
            for (int len = 1; len <= 16; ++len)
                flint_printf("\t%.1lf|%.1lf|%.1lf",
                        (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*len*100),
                        (mins_precomp[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*len*100),
                        (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*len*100));
            flint_printf("\t%.1lf|%.1lf|%.1lf",
                    (mins[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(32*32*100),
                    (mins_precomp[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(32*32*100),
                    (mins_generic[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(32*32*100));
            flint_printf("\t%.1lf|%.1lf|%.1lf",
                    (mins[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(256*256*100),
                    (mins_precomp[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(256*256*100),
                    (mins_generic[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(256*256*100));
            flint_printf("\n");
        }
        else
        {
            flint_printf("%wd", i);
            for (int len = 1; len <= 16; ++len)
                flint_printf("\t%.1lf| na |%.1lf",
                        (mins[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*len*100),
                        (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*len*100));
            flint_printf("\t%.1lf| na |%.1lf",
                    (mins[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(32*32*100),
                    (mins_generic[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(32*32*100));
            flint_printf("\t%.1lf| na |%.1lf",
                    (mins[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(256*256*100),
                    (mins_generic[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(256*256*100));
            flint_printf("\n");
        }
    }

    return 0;
}
