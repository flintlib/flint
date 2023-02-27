/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"


slong thread_pool_get_size(thread_pool_t T)
{
    slong ret;
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&T->mutex);
#endif
    ret = T->length;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&T->mutex);
#endif
    return ret;
}
