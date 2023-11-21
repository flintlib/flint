/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "acb_dft.h"
#include "profiler.h"

#define LOG 0
#define CSV 1
#define JSON 2

typedef void (*do_f) (acb_ptr y, acb_srcptr x, slong len, slong prec);

int main(int argc, char *argv[])
{
    int out;
    ulong n, nref, maxq = 5000;
    flint_rand_t state;
    slong r, nr;

    int l, nf = 5;
    do_f func[5] = {
        acb_dft_naive,
        acb_dft_crt,
        acb_dft_cyc,
        acb_dft_bluestein,
        acb_dft
    };
    char * name[6] = {
        "naive",
        "crt",
        "cyc",
        "bluestein",
        "default",
        "precomp"
    };

    int i, ni = 14;
    ulong len[14] =  { 60, 5*5*5, 131, 2*3*5*7, 256,
        3*3*3*3*3, 2*2*2*43, 503, 2*2*2*101, 3001,
        8192, 13*11*61, 10000, 101*103 };
    char * factor[14] =  { "2^2*3*5", "5^3", "131", "2*3*5*7", "2^8",
        "3^5", "2^3*43", "503", "2^3*101", "3001",
        "2^13", "11*13*61", "2^4*5^4", "101*103" };
    slong ml[14] = { 20, 10, 20, 20, 10,
        10, 10, 5, 5, 5,
        3, 3, 3, 3};

    int j, nj = 4;
    slong prec[4] = { 64, 128, 512, 1024 };
    slong mp[4] =   { 40, 20, 5, 2 };

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

    if (out == CSV)
        flint_printf("# %-12s, %7s, %7s, %7s\n","name", "prec", "len", "factor", "time");

    flint_randinit(state);

    for (j = 0; j < nj; j++)
    {

        for (i = 0; i < ni; i++)
        {
            slong k, nm;
            acb_ptr u, v;

            nm = mp[j] * ml[i];

            u = _acb_vec_init(len[i]);
            v = _acb_vec_init(len[i]);

            for (k = 0; k < len[i]; k++)
                acb_randtest_precise(u + k, state, prec[j], 4);

            if (out == LOG)
                flint_printf("%wu * DFT(%wu = %s), prec %wu....\n", nm, len[i], factor[i], prec[j]);

            for (l = 0; l <= nf; l++)
            {
                slong m;

                if (len[i] > 300 && l == 0)
                    continue;

                if (out == LOG)
                    flint_printf("%-14s ...  ", name[l]);
                else if (out == CSV)
                    flint_printf("%-12s, %7d, %8s, %7d,   ", name[l], prec[j], len[i], factor[i]);
                else if (out == JSON)
                    flint_printf("{ \"name\": \"%s\", \"prec\": %d, \"len\": %d, \"factor\": %s, \"time\": ",
                            name[l], prec[j], len[i], factor[i]);

                TIMEIT_ONCE_START

                if (l < nf)
                {
                    for (m = 0; m < nm; m++)
                        func[l](v, u, len[i], prec[j]);
                }
                else
                {
                    acb_dft_pre_t t;
                    acb_dft_precomp_init(t, len[i], prec[j]);

                    for (m = 0; m < nm; m++)
                        acb_dft_precomp(v, u, t, prec[j]);

                    acb_dft_precomp_clear(t);
                }

                TIMEIT_ONCE_STOP

                if (out == JSON)
                    flint_printf("}\n");
                else
                    flint_printf("\n");
            }

            _acb_vec_clear(u, len[i]);
            _acb_vec_clear(v, len[i]);
        }

    }
    flint_cleanup();
    return 0;
}
