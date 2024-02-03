/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"

void _nmod_vec_add_fast(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod);

#define NUMTYPES 1

void (* funcs[])(mp_ptr, mp_srcptr, mp_srcptr, slong, nmod_t) = {_nmod_vec_add /*, _nmod_vec_add_fast */};

char * str[] = {"_nmod_vec_add" /*, "_nmod_vec_add_fast" */};

typedef struct
{
    flint_bitcnt_t mod_bits;
    flint_bitcnt_t len;
    double timers[NUMTYPES];
}
info_t;

#define COUNT 10000

void sample(void * arg, ulong unused)
{
    mp_limb_t n;
    nmod_t mod;
    info_t * info = (info_t *) arg;
    flint_bitcnt_t mod_bits = info->mod_bits;
    flint_bitcnt_t len = info->len;
    slong type;
    slong ix;
    mp_ptr vec1, vec2, res;
    double * timers = info->timers;
    double start;
    FLINT_TEST_INIT(state);

    n = n_randbits(state, mod_bits);
    if (n == UWORD(0)) n++;

    nmod_init(&mod, n);

    /* warmup */
    for (type = 0; type < NUMTYPES * COUNT / 10; type++)
    {
        vec1 = _nmod_vec_init(len);
        vec2 = _nmod_vec_init(len);
        res = _nmod_vec_init(len);

        for (ix = 0; ix < len; ix++)
        {
            vec1[ix] = n_randint(state, n);
            vec2[ix] = n_randint(state, n);
        }

        for (ix = 0; ix < COUNT / len; ix++)
            funcs[type % NUMTYPES](res, vec1, vec2, len, mod);

        _nmod_vec_clear(vec1);
        _nmod_vec_clear(vec2);
    }

    /* real deal */
    for (type = 0; type < NUMTYPES * COUNT; type++)
    {
        vec1 = _nmod_vec_init(len);
        vec2 = _nmod_vec_init(len);
        res = _nmod_vec_init(len);

        for (ix = 0; ix < len; ix++)
        {
            vec1[ix] = n_randint(state, n);
            vec2[ix] = n_randint(state, n);
        }

        prof_start();
        start = clock();
        for (ix = 0; ix < COUNT / len; ix++)
            funcs[type % NUMTYPES](res, vec1, vec2, len, mod);
        timers[type % NUMTYPES] += (double) (clock() - start) / CLOCKS_PER_SEC;
        prof_stop();

        _nmod_vec_clear(vec1);
        _nmod_vec_clear(vec2);
    }

    flint_randclear(state);
}

int main(int argc, char ** argv)
{
    info_t info = {0};
    flint_bitcnt_t mod_bits;
    slong len;
    int ix;

    for (ix = 0; ix < NUMTYPES; ix++)
        printf("%25s", str[ix]);

    printf("\n");

    for (mod_bits = 3; mod_bits <= FLINT_BITS; mod_bits += 20)
    {
        info.mod_bits = mod_bits;

        for (len = 1; len <= 1200; len = 2 * len + 7)
        {
            info.len = len;

            prof_repeat(NULL, NULL, sample, (void *) &info);

            for (ix = 0; ix < NUMTYPES; ix++)
            {
                printf("%25.3lf", info.timers[ix]);
                info.timers[ix] = 0;
            }

            printf("     " WORD_FMT "d mod bits, " WORD_FMT "d len\n", mod_bits, len);
        }
    }

    return 0;
}
