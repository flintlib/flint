/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "thread_pool.h"

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
    thread_pool_handle * handles;
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

int nmod_mpoly_divides_threaded(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    slong thread_limit)
{
    slong i, * Adegs, * Bdegs;
    thread_pool_handle * handles;
    slong num_handles;
    int divides;
    TMP_INIT;

    if (B->length == 0)
    {
        flint_throw(FLINT_DIVZERO, "Exception in nmod_mpoly_divides_threaded: "
                                                   "Cannot divide by zero.\n");
    }

    if (1 != n_gcd(B->coeffs[0], ctx->ffinfo->mod.n))
    {
        flint_throw(FLINT_IMPINV, "Exception in nmod_mpoly_divides_threaded: "
                                       "Cannot invert leading coefficient.\n");
    }

    TMP_START;

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

    divides = -1;
    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS && A->length > 50)
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
            thread_pool_wake(global_thread_pool, handles[m],
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

    if (num_handles > 0)
    {
        divides = _nmod_mpoly_divides_heap_threaded(Q, A, B, ctx,
                                                         handles, num_handles);
    }
    else
    {
        divides = nmod_mpoly_divides_monagan_pearce(Q, A, B, ctx);
    }

cleanup:

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_give_back(global_thread_pool, handles[i]);
    }
    if (handles)
    {
        flint_free(handles);
    }

    TMP_END;
    return divides;
}

int nmod_mpoly_divides(
    nmod_mpoly_t Q,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_divides_threaded(Q, A, B, ctx, MPOLY_DEFAULT_THREAD_LIMIT);
}
