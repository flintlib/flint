/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#define ulong ulongxx /* fixup for Cygwin */
#include <sched.h>
#undef ulong

#include "thread_pool.h"

#if FLINT_USES_PTHREAD && defined(HAVE_PTHREAD_NP_H)
# include <pthread_np.h>
#endif

int thread_pool_restore_affinity(thread_pool_t T)
{
#if FLINT_USES_CPUSET && FLINT_USES_PTHREAD
    slong i;
    int errorno;
    thread_pool_entry_struct * D;

    /* restore affinities for workers */
    D = T->tdata;
    for (i = 0; i < T->length; i++)
    {
        errorno = pthread_setaffinity_np(D[i].pth, sizeof(cpu_set_t),
                                            (cpu_set_t *)T->original_affinity);
        if (errorno != 0)
            return errorno;
    }

    /* Restore the affinity of the calling thread. Only one affinity is
       recorded, that of the thread which called thread_pool_init, and it
       is applied here to the calling thread. This assumes that the caller
       is the thread which created the pool (for the global pool, the
       master thread, which creates it in flint_set_num_threads). If
       another thread calls this function, it receives the creator's
       affinity rather than its own original one. */
    errorno = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t),
                                        (cpu_set_t *)T->original_affinity);
    if (errorno != 0)
        return errorno;
#endif

    return 0;
}
