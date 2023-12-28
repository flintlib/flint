/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "nmod_mpoly.h"

static int _nmod_mpoly_divides_try_dense(
    slong * Adegs,
    slong * Bdegs,
    slong nvars,
    slong Alen,
    slong Blen)
{
    slong i, total_dense_size;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    total_dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, total_dense_size, total_dense_size, Adegs[i] + 1);
        if (hi != WORD(0) || total_dense_size <= WORD(0))
            return 0;
    }

    return total_dense_size < WORD(5000000)
            && total_dense_size/Alen < WORD(10);
}


typedef struct
{
    slong * degs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;
    const mpoly_ctx_struct * mctx;
    const thread_pool_handle * handles;
    slong num_handles;
}
_degrees_arg_struct;

typedef _degrees_arg_struct _degrees_arg_t[1];

static void _worker_degrees(void * varg)
{
    _degrees_arg_struct * arg = (_degrees_arg_struct *) varg;

    mpoly_degrees_si_threaded(arg->degs, arg->exps, arg->length, arg->bits,
                                    arg->mctx, arg->handles, arg->num_handles);
}

int _nmod_mpoly_divides_threaded_pool(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, * Adegs, * Bdegs;
    int divides;
    TMP_INIT;

    TMP_START;
    divides = -1;

    if (A->bits <= FLINT_BITS &&
        B->bits <= FLINT_BITS &&
        A->length > 50)
    {
        Adegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
        Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

        if (num_handles > 0)
        {
            slong m = mpoly_divide_threads(num_handles, A->length, B->length);
            _degrees_arg_t arg;

            FLINT_ASSERT(m >= 0);
            FLINT_ASSERT(m < num_handles);

            arg->degs = Bdegs;
            arg->exps = B->exps;
            arg->length = B->length;
            arg->bits = B->bits;
            arg->mctx = ctx->minfo;
            arg->handles = handles + (m + 1);
            arg->num_handles = num_handles - (m + 1);
            thread_pool_wake(global_thread_pool, handles[m], 0,
                                                         _worker_degrees, arg);
            mpoly_degrees_si_threaded(Adegs, A->exps, A->length, A->bits,
                                                   ctx->minfo, handles + 0, m);
            thread_pool_wait(global_thread_pool, handles[m]);
        }
        else
        {
            mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
            mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);
        }

        /* quick degree check */
        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            if (Adegs[i] < Bdegs[i])
            {
                nmod_mpoly_zero(Q, ctx);
                divides = 0;
                goto cleanup;
            }
        }

        if (_nmod_mpoly_divides_try_dense(Adegs, Bdegs, ctx->minfo->nvars,
                                                         A->length, B->length))
        {
            divides = nmod_mpoly_divides_dense(Q, A, B, ctx);
        }
    }

    if (divides == 0 || divides == 1)
    {
        /* have answer */
        goto cleanup;
    }

#ifdef _nmod_mpoly_divides_heap_threaded_pool
    if (num_handles > 0)
    {
        divides = _nmod_mpoly_divides_heap_threaded_pool(Q, A, B, ctx,
                                                         handles, num_handles);
    }
    else
    {
        divides = nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }
#else
    divides = nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
#endif

cleanup:

    TMP_END;
    return divides;
}

int nmod_mpoly_divides(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;
    int divides;

    if (B->length == 0)
    {
        if (A->length == 0 || nmod_mpoly_ctx_modulus(ctx) == 1)
        {
            nmod_mpoly_set(Q, A, ctx);
	        return 1;
	    }
        else
        {
    	    flint_throw(FLINT_DIVZERO, "nmod_mpoly_divides: divide by zero.");
        }
    }

    if (1 != n_gcd(B->coeffs[0], ctx->mod.n))
    {
        flint_throw(FLINT_IMPINV, "nmod_mpoly_divides: leading coefficient is not invertible.");
    }

    thread_limit = A->length/1024;

    if (A->length <= 50)
    {
        return nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

    num_handles = flint_request_threads(&handles, thread_limit);
    divides = _nmod_mpoly_divides_threaded_pool(Q, A, B, ctx, handles, num_handles);
    flint_give_back_threads(handles, num_handles);

    return divides;
}

