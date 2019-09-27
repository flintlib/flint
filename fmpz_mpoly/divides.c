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

int fmpz_mpoly_divides_threaded(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    slong thread_limit)
{
    thread_pool_handle * handles;
    slong num_handles;
    int divides;
    slong i;

    if (B->length < 2 || A->length < 2)
    {
        if (B->length == 0)
        {
            flint_throw(FLINT_DIVZERO,
                             "Divide by zero in fmpz_mpoly_divides_threaded");
        }

        if (A->length == 0)
        {
            fmpz_mpoly_zero(Q, ctx);
            return 1;
        }

        return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    handles = NULL;
    num_handles = 0;
    if (thread_limit > 1 && global_thread_pool_initialized)
    {
        slong max_num_handles;
        max_num_handles = thread_pool_get_size(global_thread_pool);
        max_num_handles = FLINT_MIN(thread_limit - 1, max_num_handles);
        if (max_num_handles > 0)
        {
            handles = (thread_pool_handle *) flint_malloc(
                                   max_num_handles*sizeof(thread_pool_handle));
            num_handles = thread_pool_request(global_thread_pool,
                                                     handles, max_num_handles);
        }
    }

    if (num_handles > 0)
    {
        divides = _fmpz_mpoly_divides_heap_threaded(Q, A, B, ctx,
                                                         handles, num_handles);
    }
    else
    {
        divides = fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_give_back(global_thread_pool, handles[i]);
    }
    if (handles)
    {
        flint_free(handles);
    }

    return divides;
}

int fmpz_mpoly_divides(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_divides_threaded(Q, A, B, ctx, MPOLY_DEFAULT_THREAD_LIMIT);
}
