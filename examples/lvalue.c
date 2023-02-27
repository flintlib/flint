/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "acb_dirichlet.h"
#include "flint/profiler.h"

int main(int argc, char *argv[])
{
    dirichlet_group_t G;
    dirichlet_char_t chi;
    acb_t s;
    acb_ptr z;
    slong prec, i, len;
    ulong q, n;
    int zfunction, deflate;

    if (argc < 2)
    {
        flint_printf("lvalue [-character q n] [-re a] [-im b] [-prec p] [-z] [-deflate] [-len l]\n\n");
        flint_printf("Print value of Dirichlet L-function at s = a+bi.\n");
        flint_printf("Default a = 0.5, b = 0, p = 53, (q, n) = (1, 0) (Riemann zeta)\n");
        flint_printf("[-z]       - compute Z(s) instead of L(s)\n");
        flint_printf("[-deflate] - remove singular term at s = 1\n");
        flint_printf("[-len l]   - compute l terms in Taylor series at s\n");
        return 1;
    }

    acb_init(s);

    zfunction = 0;
    prec = 53;
    q = 1;
    n = 0;
    len = 1;
    deflate = 0;
    acb_set_d(s, 0.5);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-prec"))
        {
            prec = atol(argv[i+1]);
        }
    }

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-character"))
        {
            q = atol(argv[i+1]);
            n = atol(argv[i+2]);
        }
        else if (!strcmp(argv[i], "-len"))
        {
            len = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-z"))
        {
            zfunction = 1;
        }
        else if (!strcmp(argv[i], "-deflate"))
        {
            deflate = 1;
        }
        else if (!strcmp(argv[i], "-re"))
        {
            arb_set_str(acb_realref(s), argv[i+1], prec);
        }
        else if (!strcmp(argv[i], "-im"))
        {
            arb_set_str(acb_imagref(s), argv[i+1], prec);
        }
    }

    if (n_gcd(q, n) != 1)
    {
        flint_printf("need gcd(q,n) = 1 to define a character\n");
        flint_abort();
    }

    z = _acb_vec_init(len);

    dirichlet_group_init(G, q);
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, n);

    TIMEIT_ONCE_START
    if (zfunction)
        acb_dirichlet_hardy_z(z, s, G, chi, len, prec);
    else
        acb_dirichlet_l_jet(z, s, G, chi, deflate, len, prec);

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            flint_printf("%s(s) = ", zfunction ? "Z" : "L");
        else if (i == 1)
            flint_printf("%s'(s) = ", zfunction ? "Z" : "L");
        else
            flint_printf("[x^%wd] %s(s+x) = ", i, zfunction ? "Z" : "L");
        acb_printn(z + i, prec * 0.333 + 1, 0);
        flint_printf("\n");
    }

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    _acb_vec_clear(z, len);
    acb_clear(s);
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);
    flint_cleanup();
    return EXIT_SUCCESS;
}

