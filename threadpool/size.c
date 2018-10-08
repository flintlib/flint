/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"


slong threadpool_size(threadpool_t T)
{
    slong ret;
    pthread_mutex_lock(&T->mutex);
    ret = T->length;
    pthread_mutex_unlock(&T->mutex);
    return ret;
}
