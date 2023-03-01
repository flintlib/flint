/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "acb_dirichlet.h"
#include "profiler.h"

#define LOG 0
#define CSV 1
#define JSON 2

typedef void (*do_f) (acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

int main(int argc, char *argv[])
{
    int out;
    ulong n, nref, maxq = 5000;
    ulong * rand;
    flint_rand_t state;
    slong r, nr;

    int l, nf = 5;
    do_f func[5] = { acb_dirichlet_gauss_sum_naive,
        acb_dirichlet_gauss_sum_factor,
        acb_dirichlet_gauss_sum_theta,
        acb_dirichlet_gauss_sum,
        acb_dirichlet_gauss_sum };
    char * name[5] = { "naive", "factor", "theta", "default", "precomp" };

    int i, ni = 5;
    ulong qmin[5] =  {  3,  50, 500, 1000, 10000 };
    ulong qmax[5] =  { 50, 100, 599, 1050, 10030 };

    int j, nj = 3;
    slong prec[3] = { 64, 128, 1024 };

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
        flint_printf("# %-12s, %7s, %7s, %7s, %7s\n","name", "prec", "qmin", "qmax", "time");

    for (j = 0; j < nj; j++)
    {

    for (i = 0; i < ni; i++)
    {

        if (out == LOG)
            flint_printf("G_q(a) at prec %wu for all %wu <= q <= %wu....\n", prec[j], qmin[i], qmax[i]);

        for (l = 0; l < nf; l++)
        {
            ulong q;

            if (qmin[i] > 2000 && l == 0)
                continue;

            if (out == LOG)
                flint_printf("%-14s ...  ", name[l]);
            else if (out == CSV)
                flint_printf("%-12s, %7d, %7d, %7d,   ", name[l], prec[j], qmin[i], qmax[i]);
            else if (out == JSON)
                flint_printf("{ \"name\": \"%s\", \"prec\": %d, \"qmin\": %d, \"qmax\": %d, \"time\": ",
                        name[l], prec[j], qmin[i], qmax[i]);

            TIMEIT_ONCE_START

                for (q = qmin[i]; q <= qmax[i]; q++)
                {
                    dirichlet_group_t G;
                    dirichlet_char_t chi;
                    acb_t res;

                    if (q % 4 == 2)
                        continue;

                    dirichlet_group_init(G, q);
                    if (l == 4)
                        dirichlet_group_dlog_precompute(G, 1);
                    dirichlet_char_init(chi, G);
                    dirichlet_char_first_primitive(chi, G);
                    acb_init(res);

                        do {
                            func[l](res, G, chi, prec[j]);
                        } while (dirichlet_char_next_primitive(chi, G) >= 0);

                    acb_clear(res);
                    dirichlet_char_clear(chi);
                    if (l == 4)
                        dirichlet_group_dlog_clear(G);
                    dirichlet_group_clear(G);
                }

            TIMEIT_ONCE_STOP

                if (out == JSON)
                    flint_printf("}\n");
                else
                    flint_printf("\n");
        }

    }

    }
    flint_cleanup();
    return EXIT_SUCCESS;
}
