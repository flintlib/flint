/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "dlog.h"
#include "profiler.h"

#define NUMPRIMES 400
#define LOG 0
#define CSV 1
#define JSON 2

typedef void (*log_f) (ulong p, ulong a, ulong num);

void
flog_table(ulong p, ulong a, ulong num)
{
    int k;
    dlog_table_t t;
    dlog_table_init(t, a, p);
    for (k = 1; k < num; k++)
    {
        if (k % p == 0) continue;
        dlog_table(t, k % p);
    }
    dlog_table_clear(t);
}
void
flog_bsgs(ulong p, ulong a, ulong num)
{
    int k;
    dlog_bsgs_t t;
    dlog_bsgs_init(t, a, p, p-1, dlog_bsgs_size(p, num));
    for (k = 1; k < num; k++)
    {
        if (k % p == 0) continue;
        dlog_bsgs(t, k % p);
    }
    dlog_bsgs_clear(t);
}
void
flog_rho(ulong p, ulong a, ulong num)
{
    int k;
    dlog_rho_t t;
    dlog_rho_init(t, a, p, p-1);
    for (k = 1; k < num; k++)
    {
        if (k % p == 0) continue;
        dlog_rho(t, k % p);
    }
    dlog_rho_clear(t);
}
void
flog_crt(ulong p, ulong a, ulong num)
{
    int k;
    dlog_crt_t t;
    dlog_crt_init(t, a, p, p-1, num);
    for (k = 1; k < num; k++)
    {
        if (k % p == 0) continue;
        dlog_crt(t, k % p);
    }
    dlog_crt_clear(t);
}
void
flog_gen(ulong p, ulong a, ulong num)
{
    int k;
    dlog_precomp_t t;
    dlog_precomp_n_init(t, a, p, p-1, num);
    for (k = 1; k < num; k++)
    {
        if (k % p == 0) continue;
        dlog_precomp(t, k % p);
    }
    dlog_precomp_clear(t);
}

int main(int argc, char *argv[])
{
    int out = LOG;
    slong iter, k, nv, nref, r, nr;
    ulong minq, maxq;
    ulong * rand;
    int nbits, nl = 5;
    int l[5] = { 1, 10, 100, 1000 , 5000};

    int nf = 4;
    log_f func[4] = { flog_table, flog_bsgs, flog_crt, flog_gen };
    char * n[4] = { "table", "bsgs", "crt", "generic" };

    int np = NUMPRIMES;

    flint_rand_t state;

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

    flint_randinit(state);
    for (nbits = 10; nbits <= 40; nbits += 5)
    {

        int i;
        ulong p[NUMPRIMES], a[NUMPRIMES];

        if (nbits >= 25)
            np /= 2;

        for (i=0; i < np; i++)
        {
            p[i] = n_randprime(state, nbits, 0);
            a[i] = n_primitive_root_prime(p[i]);
        }

        for (i = 0; i < nl; i++)
        {
            int f;

            if (out == LOG)
                flint_printf("%d * logs mod primes of size %d.....\n", l[i], nbits);

            for (f = 0; f < nf; f++)
            {
                int j;

                /* skip useless */
                if (f == 0 && nbits >= 20)
                    continue;
                if (f == 1 && nbits >= 30 && l[i] > 10)
                    continue;
                if (out == LOG)
                {
                    flint_printf("%-20s...   ",n[f]);
                    fflush(stdout);
                }
                else if (out == CSV)
                    flint_printf("%-8s, %2d, %4d, %3d, ",n[f],nbits,l[i],np);
                else if (out == JSON)
                    flint_printf("{ \"name\": \"%s\", \"bits\": %d, \"nlogs\": %d, \"nprimes\": %d, \"time\": ",
                            n[f],nbits,l[i],np);

                TIMEIT_ONCE_START
                    for (j = 0; j < np; j ++)
                        (func[f])(p[j], a[j], l[i]);
                TIMEIT_ONCE_STOP

                    if (out == JSON)
                        flint_printf("}\n");
                    else
                        flint_printf("\n");
            }
        }
    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
