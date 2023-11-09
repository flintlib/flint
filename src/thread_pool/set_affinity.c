/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <sched.h>

#include "thread_pool.h"

#if FLINT_USES_PTHREAD && defined(HAVE_PTHREAD_NP_H)
# include <pthread_np.h>
#endif

/*
    cpus[0] is cpu number for main thread
    cpus[1], ..., cpus[length - 1] are cpu numbers for elements of pool
*/
int thread_pool_set_affinity(thread_pool_t T, int * cpus, slong length)
{
#if FLINT_USES_CPUSET && FLINT_USES_PTHREAD
    slong i;
    int errorno;
    cpu_set_t mask;
    thread_pool_entry_struct * D;

    if (length <= 0)
        return 0;

    /* set affinities for workers */
    D = T->tdata;
    for (i = 0; i + 1 < length && i < T->length; i++)
    {
        CPU_ZERO(&mask);
        CPU_SET(cpus[i + 1] % CPU_SETSIZE, &mask);
        errorno = pthread_setaffinity_np(D[i].pth, sizeof(cpu_set_t), &mask);
        if (errorno != 0)
            return errorno;
    }

    /* set affinity for main thread */
    CPU_ZERO(&mask);
    CPU_SET(cpus[0] % CPU_SETSIZE, &mask);
    errorno = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask);
    if (errorno != 0)
        return errorno;
#endif

    return 0;
}
