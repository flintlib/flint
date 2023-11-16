/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "dirichlet.h"
#include "profiler.h"

#define LOG 0
#define CSV 1
#define JSON 2

typedef ulong (*do_f) (ulong q1, ulong q2);

static ulong
do_gcd(ulong q1, ulong q2)
{
    ulong n, q, k;

    for (n = 0, q = q1; q <= q2; q++)
        for (k = 1; k < q; k++)
            n += (n_gcd(k, q) == 1);

    return n;
}

static ulong
do_char(ulong q1, ulong q2)
{
    ulong n, q;

    for (n = 0, q = q1; q <= q2; q++)
    {
        dirichlet_group_t G;
        dirichlet_char_t x;

        dirichlet_group_init(G, q);
        dirichlet_char_init(x, G);

        dirichlet_char_one(x, G);
        n++;
        for (; dirichlet_char_next(x, G) >= 0; n++);

        dirichlet_char_clear(x);
        dirichlet_group_clear(G);
    }

    return n;
}

static ulong
do_gcdpluscond(ulong q1, ulong q2)
{
    ulong n, q, k;

    for (n = 0, q = q1; q <= q2; q++)
    {
        dirichlet_group_t G;
        dirichlet_group_init(G, q);

        for (k = 1; k < q; k++)
        {
            /* known factors -> faster gcd */
            slong i;
            if (G->q_even > 1 && k % 2 == 0)
                continue;
            for (i = G->neven; i < G->num; i++)
                if (k % G->P[i].p == 0)
                    break;

            if (i == G->num)
                dirichlet_conductor_ui(G, k);
        }

        n += G->phi_q;
        dirichlet_group_clear(G);

    }

    return n;
}

static ulong
do_charpluscond(ulong q1, ulong q2)
{
    ulong n, q, k;

    for (n = 0, q = q1; q <= q2; q++)
    {
        dirichlet_group_t G;
        dirichlet_char_t x;

        dirichlet_group_init(G, q);
        dirichlet_char_init(x, G);

        dirichlet_char_one(x, G);
        n++;

        for (; dirichlet_char_next(x, G) >= 0; n++)
            dirichlet_conductor_char(G, x);

        dirichlet_char_clear(x);
        dirichlet_group_clear(G);

    }

    return n;
}

static ulong
do_charplusorder(ulong q1, ulong q2)
{
    ulong n, q, k;

    for (n = 0, q = q1; q <= q2; q++)
    {
        dirichlet_group_t G;
        dirichlet_char_t x;

        dirichlet_group_init(G, q);
        dirichlet_char_init(x, G);

        dirichlet_char_one(x, G);
        n++;

        for (; dirichlet_char_next(x, G) >= 0; n++)
            dirichlet_order_char(G, x);

        dirichlet_char_clear(x);
        dirichlet_group_clear(G);

    }

    return n;
}

int main(int argc, char *argv[])
{
    int out;
    ulong n, nref, maxq = 5000;

    int l, nf = 5;
    do_f func[5] = {
        do_gcd,
        do_char,
        do_gcdpluscond,
        do_charpluscond,
        do_charplusorder
    };
    char * name[5] = {
        "gcd",
        "char",
        "gcd + cond",
        "char + cond",
        "char + order"
    };

    int i, ni = 5;
    ulong qmin[5] = { 2,   1000, 10000, 100000, 1000000 };
    ulong qmax[5] = { 500, 3000, 11000, 100100, 1000010 };

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
        flint_printf("# %-12s, %7s, %7s, %7s\n","name", "qmin", "qmax", "time");

    for (i = 0; i < ni; i++)
    {

        if (out == LOG)
        {
            flint_printf("loop over all (Z/q)* for %wu<=q<=%wu\n", qmin[i], qmax[i]);
        }

        for (l = 0; l < nf; l++)
        {

            if (out == LOG)
                flint_printf("%-14s ...  ",name[l]);
            else if (out == CSV)
                flint_printf("%-12s, %7d, %7d,   ",name[l],qmin[i],qmax[i]);
            else if (out == JSON)
                flint_printf("{ \"name\": \"%s\", \"qmin\": %d, \"qmax\": %d, \"time\": ",
                        name[l],qmin[i],qmax[i]);

            TIMEIT_ONCE_START
                (func[l])(qmin[i], qmax[i]);
            TIMEIT_ONCE_STOP

            if (l == 0)
                nref = n;
            else if (n != nref)
            {
                flint_printf("FAIL: wrong number of elements %wu != %wu\n\n",n, nref);
                flint_abort();
            }

            if (out == JSON)
                flint_printf("}\n");
            else
                flint_printf("\n");
        }

    }

    flint_cleanup();
    return 0;
}
