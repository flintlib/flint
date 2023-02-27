/* This file is public domain. Author: Pascal Molin. */

#include <string.h>
#include <math.h>
#include "flint/profiler.h"
#include "acb_dirichlet.h"

static int usage(char *argv[])
{
    printf("Computes central values (s = 0.5) of Dirichlet L-functions.\n\n");
    printf("usage: %s [--quiet] [--check] [--prec <bits>] qmin qmax\n", argv[0]);
    return 1;
}

int main(int argc, char *argv[])
{
    int i, check = 0, out = 1;
    slong prec = 100, digits = 30;
    ulong qmin, qmax, q;
    acb_t s;
    acb_dirichlet_hurwitz_precomp_t pre;

    if (argc < 3)
        return usage(argv);

    for (i = 1; i < argc - 2; i++)
    {
        if (!strcmp(argv[i], "--quiet"))
            out = 0;
        else if (!strcmp(argv[i], "--check"))
            check = 1;
        else if (!strcmp(argv[i], "--prec"))
        {
            i++;
            prec = atol(argv[i]);
            digits = floor(prec * 0.3);
        }
    }

    if (argc < i + 2)
        return usage(argv);

    qmin = atol(argv[i]);
    qmax = atol(argv[i+1]);

    fflush(stdout);

    acb_init(s);

    TIMEIT_ONCE_START

    acb_one(s);
    acb_div_si(s, s, 2, prec);

    acb_dirichlet_hurwitz_precomp_init_num(pre, s, 0,
        (qmax - qmin + 1) * 0.5 * qmax, prec);

    for (q = qmin; q <= qmax; q++)
    {
        ulong k;
        dirichlet_group_t G;
        dirichlet_char_t x;
        acb_ptr z;

        if (q % 4 == 2)
            continue;

        dirichlet_group_init(G, q);
        dirichlet_char_init(x, G);

        z =  _acb_vec_init(G->phi_q);

        acb_dirichlet_l_vec_hurwitz(z, s, pre, G, prec);

        if (out || check)
        {
            k = 0;
            dirichlet_char_one(x, G);
            while (dirichlet_char_next(x, G) >= 0)
            {
                k++;
                if (dirichlet_conductor_char(G,x) < q)
                    continue;
                if (acb_contains_zero(z + k))
                {
                    flint_printf("Value could be zero!\n");
                    flint_printf("%wu,%wu: ", q, x->n);
                    acb_printn(z + k, digits, 0);
                    flint_printf("\n");
                    flint_abort();
                }

                if (out)
                {
                    flint_printf("%wu,%wu: ", q, x->n);
                    acb_printn(z + k, digits, 0);
                    flint_printf("\n");
                }
            }
        }

        _acb_vec_clear(z, G->phi_q);
        dirichlet_char_clear(x);
        dirichlet_group_clear(G);
    }

    acb_dirichlet_hurwitz_precomp_clear(pre);

    acb_clear(s);

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}

