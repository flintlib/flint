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

typedef void (*do_f) (ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv);

static void
do_empty(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
    return;
}

static void
do_dlog_primeloop(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
    slong k, l;
	nmod_t order;
	nmod_init(&order, dirichlet_order_char(G, chi));

    for (k = 0; k < nv; k++)
        v[k] = 0;

    for (l = G->neven; l < G->num; l++)
    {
        dirichlet_prime_group_struct P = G->P[l];
        dlog_vec_loop_add(v, nv, P.g, chi->log[l], P.pe, P.phi.n, order);
    }
    dirichlet_vec_set_null(v, G, nv);
}

static void
do_eratos(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
	slong k, p, pmax;
	nmod_t order;
	n_primes_t iter;

	n_primes_init(iter);
	nmod_init(&order, dirichlet_order_char(G, chi));

	pmax = (nv < G->q) ? nv : G->q;
	v[1] = 0;

	while ((p = n_primes_next(iter)) < pmax)
	{
		if (G->q % p == 0)
		{
			for (k = p; k < nv; k += p)
				v[k] = DIRICHLET_CHI_NULL;
		}
		else
		{
			long chip;
			chip = dirichlet_chi(G, chi, p);

			for (k = p; k < nv; k += p)
				if (v[k] != -1)
					 v[k] = nmod_add(v[k], chip, order);
		}
	}

	n_primes_clear(iter);
}

int main(int argc, char *argv[])
{
    int out;
    ulong n, nref, maxq = 5000;
    ulong * rand;
    flint_rand_t state;
    slong r, nr;

    int l, nf = 9;
    do_f func[9] = {
        do_empty,
        dirichlet_chi_vec_loop,
        do_dlog_primeloop,
        dirichlet_chi_vec_primeloop,
        do_eratos,
        dirichlet_chi_vec,
        dirichlet_chi_vec,
        dirichlet_chi_vec,
        dirichlet_chi_vec };
    char * name[9] = {
        "char only",
        "big loop",
        "prime loops",
        "prime dlog_vec",
        "manual eratos",
        "default",
        "precomp 1",
        "precomp 20",
        "precomp 100" };

    int i, ni = 5;
    ulong qmin[5] =  {   2, 1000, 3000, 10000, 100000 };
    ulong qmax[5] =  { 500, 2000, 5000, 12000, 100500 };

    int j, nj = 3;
    slong nv[3] = { 50, 300, 2000 };

    nr = 20;

    flint_randinit(state);

    rand = flint_malloc(nr * sizeof(ulong));

    for (r = 0; r < nr; r++)
        rand[r] = n_randprime(state, 42, 0);

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
        flint_printf("# %-12s, %7s, %7s, %7s, %7s\n","name", "num", "qmin", "qmax", "time");

    for (j = 0; j < nj; j++)
    {

        ulong * v;
        v = flint_malloc(nv[j] * sizeof(ulong));

        for (i = 0; i < ni; i++)
        {

            if (out == LOG)
                flint_printf("%wu * ui_chi(rand, 1..%wu) for all %wu <= q <= %wu....\n", nr, nv[j], qmin[i], qmax[i]);

            for (l = 0; l < nf; l++)
            {
                ulong q;

                /* eratos too slow */
                if (l == 4 && i > 2)
                    continue;

                if (out == LOG)
                    flint_printf("%-14s ...  ", name[l]);
                else if (out == CSV)
                    flint_printf("%-12s, %7d, %7d, %7d,   ", name[l], nv[j], qmin[i], qmax[i]);
                else if (out == JSON)
                    flint_printf("{ \"name\": \"%s\", \"num\": %d, \"qmin\": %d, \"qmax\": %d, \"time\": ",
                            name[l], nv[j], qmin[i], qmax[i]);

                TIMEIT_ONCE_START

                for (q = qmin[i]; q <= qmax[i]; q++)
                {
                    dirichlet_group_t G;
                    dirichlet_char_t chi;

                    dirichlet_group_init(G, q);
                    dirichlet_char_init(chi, G);

                    if (l >= 6)
                        dirichlet_group_dlog_precompute(G, (l == 6) ? 1 : (l==7) ? 20 : 100);

                    for (r = 0; r < nr; r++)
                    {
                        dirichlet_char_log(chi, G, rand[r] % q);
                        func[l](v, G, chi, nv[j]);
                    }

                    if (l >= 6)
                        dirichlet_group_dlog_clear(G);

                    dirichlet_char_clear(chi);
                    dirichlet_group_clear(G);
                }

                TIMEIT_ONCE_STOP

                    if (out == JSON)
                        flint_printf("}\n");
                    else
                        flint_printf("\n");
            }

        }
        flint_free(v);

    }

    flint_free(rand);
    flint_randclear(state);

    flint_cleanup();
    return 0;
}
