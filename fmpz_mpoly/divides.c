/*
    Copyright (C) 2018-2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "thread_pool.h"

int fmpz_mpoly_divides(fmpz_mpoly_t Q, const fmpz_mpoly_t A,
                              const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    int ret;

    if (!global_thread_pool_initialized)
    {
        ret = fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }
    else
    {
        slong max_num_workers = thread_pool_get_size(global_thread_pool);
        if (A->length > 64*max_num_workers)
        {
            ret = fmpz_mpoly_divides_heap_threaded(Q, A, B, ctx);
        }
        else
        {
            ret = fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
        }
    }

    return ret;
}
