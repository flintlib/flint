/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"

slong flint_thread_pool_num_available(thread_pool_t T)
{
    slong i, num = 0;
    thread_pool_entry_struct * D;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&T->mutex);
#endif

    D = T->tdata;
    if (T->length > 0)
    {
        for (i = 0; i < T->length; i++)
        {
            num += D[i].available;
        }
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&T->mutex);
#endif

    return num;
}

slong flint_get_num_available_threads(void)
{
    if (global_thread_pool_initialized)
        return flint_thread_pool_num_available(global_thread_pool) + 1;

    return flint_get_num_threads();
}
