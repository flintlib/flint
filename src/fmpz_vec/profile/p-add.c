/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

#define NUMTYPES 2

typedef struct
{
    flint_bitcnt_t bits;
    flint_bitcnt_t len;
    ulong timers[NUMTYPES];
}
info_t;

#define COUNT 100000

void profiler(void * arg)
{
    info_t * info = (info_t *) arg;
    flint_bitcnt_t bits = info->bits;
    flint_bitcnt_t len = info->len;
    slong ix;
    fmpz * r1, * r2, * vec;
    ulong * timers = info->timers;
    ulong start;
    FLINT_TEST_INIT(state);

    /* warmup */
    for (ix = 0; ix < COUNT / (8 * len); ix++)
    {
        r1 = _fmpz_vec_init(len);
        r2 = _fmpz_vec_init(len);
        vec = _fmpz_vec_init(len);

        _fmpz_vec_randtest(r1, state, len, bits);
        _fmpz_vec_randtest(vec, state, len, bits);
        _fmpz_vec_set(r2, r1, len);

        _fmpz_vec_add(r1, r1, vec, len);
        _fmpz_vec_inplace_add(r2, vec, len);

        _fmpz_vec_clear(r1, len);
        _fmpz_vec_clear(r2, len);
        _fmpz_vec_clear(vec, len);
    }

    /* real deal */
    for (ix = 0; ix < COUNT / len; ix++)
    {
        r1 = _fmpz_vec_init(len);
        r2 = _fmpz_vec_init(len);
        vec = _fmpz_vec_init(len);

        _fmpz_vec_randtest(r1, state, len, bits);
        _fmpz_vec_randtest(vec, state, len, bits);
        _fmpz_vec_set(r2, r1, len);

        start = clock();
        _fmpz_vec_add(r1, r1, vec, len);
        timers[0] += clock() - start;

        start = clock();
        _fmpz_vec_inplace_add(r2, vec, len);
        timers[1] += clock() - start;

        _fmpz_vec_clear(r1, len);
        _fmpz_vec_clear(r2, len);
        _fmpz_vec_clear(vec, len);
    }

    flint_randclear(state);
}

int main(int argc, char ** argv)
{
    info_t info = {0};
    flint_bitcnt_t bits;
    slong len;
    int ix;

    printf("%15s%15s\n", "_fmpz_vec_add", "_fmpz_vec_inplace_add");

    for (bits = 3; bits <= 3 * FLINT_BITS; bits += 20)
    {
        info.bits = bits;

        for (len = 1; len <= 8000; len = 2 * len + 1)
        {
            info.len = len;

            profiler(&info);

            for (ix = 0; ix < NUMTYPES; ix++)
            {
                printf("%15.3lf", 1000.0 * (double) info.timers[ix] / CLOCKS_PER_SEC);
                info.timers[ix] = 0;
            }

            printf("   " WORD_FMT "d bits, " WORD_FMT "d len\n", bits, len);
        }
    }

    return 0;
}
