/*
   Generate irreducible polynomial of minimal weight over GF(p)
   for degrees nmin <= n <= nmax.

   This file is public domain. Author: Fredrik Johansson.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "flint/thread_support.h"
#include "flint/ulong_extras.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_poly_factor.h"
#include "flint/profiler.h"

typedef struct
{
    nmod_poly_struct * f;
    ulong nstart;
}
workinfo_t;

void
worker(slong i, void * work)
{
    nmod_poly_struct * f = ((workinfo_t *) work)->f;
    ulong nstart = ((workinfo_t *) work)->nstart;
    nmod_poly_minimal_irreducible(f + i, nstart + i);
}

int
main(int argc, char * argv[])
{
    slong i;
    int num_threads = 1;
    ulong p, n, na, nb, nmin, nmax;
    slong batch_size;
    workinfo_t work;

    if (argc < 4)
    {
        flint_printf("usage: minimal_irreducibles [-threads t] p nmin nmax\n");
        return 1;
    }

    p = 0;
    nmin = 0;
    nmax = 0;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atoi(argv[i+1]);
            flint_set_num_threads(num_threads);
            i++;
        }
        else if (i == argc - 3)
            p = atol(argv[i]);
        else if (i == argc - 2)
            nmin = atol(argv[i]);
        else if (i == argc - 1)
            nmax = atol(argv[i]);
    }

    flint_printf("p = %wu, nmin = %wu, nmax = %wu\n", p, nmin, nmax);

    if (!n_is_prime(p))
        flint_throw(FLINT_ERROR, "p = %wu is not prime\n", p);

    if (num_threads == 1)
        batch_size = 1;
    else
        batch_size = FLINT_MIN(10 * num_threads, nmax - nmin + 1);
    work.f = flint_malloc(sizeof(nmod_poly_struct) * batch_size);
    for (i = 0; i < batch_size; i++)
        nmod_poly_init(work.f + i, p);

    TIMEIT_ONCE_START;
    na = nmin;
    while (1)
    {
        work.nstart = na;
        nb = FLINT_MIN(na + batch_size - 1, nmax);
        flint_parallel_do(worker, &work, nb - na + 1, 0, 0);
        for (n = na; n <= nb; n++)
            flint_printf("%{nmod_poly}\n", work.f + n - na);
        na = nb + 1;
        if (na > nmax)
            break;
    }
    TIMEIT_ONCE_STOP;
    print_memory_usage();

    for (i = 0; i < batch_size; i++)
        nmod_poly_clear(work.f + i);
    flint_free(work.f);

    flint_cleanup_master();
    return 0;
}

