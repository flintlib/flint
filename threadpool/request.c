/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"


slong threadpool_request(threadpool_t T, threadpool_threadhandle * out,
                                                               slong requested)
{
    slong i, ret = 0;
    tpentry_struct * D;

    if (requested <= 0)
        return 0;

    pthread_mutex_lock(&T->mutex);

    D = T->tdata;
    if (T->length > 0)
    {
        for (i = 0; i < T->length; i++)
        {
            if (D[i].available == 1)
            {
                D[i].available = 0;
                out[ret] = i;
                ret++;
                if (ret >= requested)
                    break;
            }
        }
    }

    pthread_mutex_unlock(&T->mutex);

    return ret;
}
