/*
    Copyright (C) 2018-2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

#ifdef _fmpz_mpoly_divides_heap_threaded_pool

#include "thread_support.h"

int fmpz_mpoly_divides(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    thread_pool_handle * handles;
    slong num_handles;
    int divides;
    slong thread_limit;

    thread_limit = A->length/1024;

    if (B->length < 2 || A->length < 2)
    {
        if (B->length == 0)
        {
            flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides");
        }

        if (A->length == 0)
        {
            fmpz_mpoly_zero(Q, ctx);
            return 1;
        }

        return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    num_handles = flint_request_threads(&handles, thread_limit);

    divides = (num_handles > 0)
            ? _fmpz_mpoly_divides_heap_threaded_pool(Q, A, B, ctx,
                                                          handles, num_handles)
            : fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);

    flint_give_back_threads(handles, num_handles);

    return divides;
}
#else
int fmpz_mpoly_divides(
    fmpz_mpoly_t Q,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides");
    }

    if (A->length == 0)
    {
        fmpz_mpoly_zero(Q, ctx);
        return 1;
    }

    return fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx);
}
#endif
