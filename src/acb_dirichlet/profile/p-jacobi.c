/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "acb_dirichlet.h"
#include "profiler.h"

#define LOG 0
#define CSV 1
#define JSON 2

typedef void (*do_f) (acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);

int main(int argc, char *argv[])
{
    int out;
    ulong n, nref, maxq = 5000;
    ulong * rand;
    flint_rand_t state;
    slong r, nr;

    int l, nf = 4;
    do_f func[4] = { acb_dirichlet_jacobi_sum_naive,
        acb_dirichlet_jacobi_sum_factor,
        acb_dirichlet_jacobi_sum_gauss,
        acb_dirichlet_jacobi_sum };
    char * name[4] = { "naive", "factor", "gauss", "default" };

    int i, ni = 4;
    ulong qmin[4] =  {  3,  51, 501, 10001 };
    ulong qmax[4] =  { 30, 100, 550, 10050 };

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
            flint_printf("J_q(a,b) at prec %wu for all %wu <= q <= %wu....\n", prec[j], qmin[i], qmax[i]);

        for (l = 0; l < nf; l++)
        {
            ulong q;

            if (out == LOG)
                flint_printf("%-14s ...  ", name[l]);
            else if (out == CSV)
                flint_printf("%-12s, %7d, %7d, %7d,   ", name[l], prec[j], qmin[i], qmax[i]);
            else if (out == JSON)
                flint_printf("{ \"name\": \"%s\", \"prec\": %d, \"qmin\": %d, \"qmax\": %d, \"time\": ",
                        name[l], prec[j], qmin[i], qmax[i]);

            TIMEIT_ONCE_START

                for (q = qmin[i]; q <= qmax[i]; q+=2)
                {
                    dirichlet_group_t G;
                    dirichlet_char_t chi1, chi2;
                    acb_t res;

                    dirichlet_group_init(G, q);
                    dirichlet_char_init(chi1, G);
                    dirichlet_char_init(chi2, G);
                    acb_init(res);

                    do {
                        slong c = 0;

                        dirichlet_char_set(chi2, G, chi1);
                        do {
                            func[l](res, G, chi1, chi2, prec[j]);
                        } while (dirichlet_char_next(chi2, G) >= 0 && c++ < 30);

                    } while (dirichlet_char_next(chi1, G) >= 0);

                    acb_clear(res);
                    dirichlet_char_clear(chi1);
                    dirichlet_char_clear(chi2);
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
    return 0;
}
