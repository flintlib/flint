/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "nmod.h"
#include "nmod_vec.h"
#include "dlog.h"
#include "profiler.h"

#define NPRIMES 640
#define LOG 0
#define CSV 1
#define JSON 2

typedef void (*vec_f) (ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order);

void
f_empty(ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order)
{
    return;
}

int main(int argc, char *argv[])
{
    int out = LOG;
    int i, ni = 8;
    int bits[9] = { 10, 15, 20, 25, 30, 35, 40, 45, 50 };

    int j, nj = 6;
    ulong * v;
    ulong nv[6] = { 50, 200, 1000, 2000, 10000, 30000 };

    int k, np = NPRIMES;
    nmod_t * p;
    ulong * a;

    int l, nf = 5;
    vec_f func[5] = { f_empty, dlog_vec_loop, dlog_vec_eratos, dlog_vec_sieve, dlog_vec };
    char * n[5] = { "empty", "loop", "eratos", "sieve", "default" };

    flint_rand_t state;
    nmod_t order;

    nmod_init(&order, 100);
    p = flint_malloc(np * sizeof(nmod_t));
    a = flint_malloc(np * sizeof(ulong));

    flint_randinit(state);

    if (argc < 2)
        out = LOG;
    else if (!strcmp(argv[1], "json"))
        out = JSON;
    else if (!strcmp(argv[1], "csv"))
        out = CSV;
    else if (!strcmp(argv[1], "log"))
        out = LOG;
    else
    {
        printf("usage: %s [log|csv|json]\n", argv[0]);
        flint_abort();
    }

    for (i = 0; i < ni; i++)
    {

        for (k = 0; k < np; k++)
        {
            nmod_init(&p[k], n_randprime(state, bits[i], 0));
            a[k] = n_primitive_root_prime(p[k].n);
        }

        for (j = 0; j < nj; j++)
        {

            v = flint_malloc(nv[j] * sizeof(ulong));

            if (out == LOG)
            {
                flint_printf("log(1..%wu) mod %d primes of size %d bits....\n", nv[j], np, bits[i]);
                fflush(stdout);
            }

            for (l = 0; l < nf; l++)
            {
                if (l == 1 && i > 2)
                    continue;
                if (l == 2 && i > 5)
                    continue;

                if (out == LOG)
                {
                    flint_printf("%-20s...   ",n[l]);
                    fflush(stdout);
                }
                else if (out == CSV)
                {
                    flint_printf("%-8s, %2d, %4d, %3d, ",n[l],bits[i],nv[j],np);
                }
                else if (out == JSON)
                {
                    flint_printf("{ \"name\": \"%s\", \"bits\": %d, \"nv\": %d, \"nprimes\": %d, \"time\": ",
                            n[l],bits[i],nv[j],np);
                }

                TIMEIT_ONCE_START

                for (k = 0; k < np; k++)
                {
                    int kk;

                    for (kk=0; kk < nv[j]; kk++)
                        v[kk] = 0;

                    (func[l])(v, nv[j], a[k], 1, p[k], p[k].n - 1, order);
                }

                TIMEIT_ONCE_STOP

                if (out == JSON)
                    flint_printf("}\n");
                else
                    flint_printf("\n");
            }
            flint_free(v);
        }
        np /= 2;
    }

    flint_free(p);
    flint_free(a);
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
