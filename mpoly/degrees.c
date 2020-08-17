/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

void mpoly_degrees_si(
    slong * user_degs,
    const ulong * poly_exps,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * max_fields;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            user_degs[i] = -WORD(1);
        return;
    }

    TMP_START;

    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields + i);

    mpoly_max_fields_fmpz(max_fields, poly_exps, len, bits, mctx);

    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)user_degs, max_fields, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields + i);

    TMP_END;
}


typedef struct
{
    slong * degs;
    const ulong * start;
    slong length;
    flint_bitcnt_t bits;
    const mpoly_ctx_struct * mctx;
}
_degrees_si_arg_struct;

static void _worker_degrees_si(void * varg)
{
    _degrees_si_arg_struct * arg = (_degrees_si_arg_struct *) varg;
    mpoly_degrees_si(arg->degs, arg->start, arg->length, arg->bits, arg->mctx);
}

void mpoly_degrees_si_threaded(
    slong * user_degs,
    const ulong * poly_exps,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j;
    slong num_threads;
    _degrees_si_arg_struct * args;
    slong start, stop;
    slong N = mpoly_words_per_exp(bits, mctx);
    slong * degs_array;

    if (len == 0)
    {
        for (j = 0; j < mctx->nvars; j++)
            user_degs[j] = -WORD(1);
        return;
    }

    num_threads = num_handles + 1;

    degs_array = (slong *) flint_malloc(num_threads*mctx->nvars*sizeof(slong));
    args = (_degrees_si_arg_struct *) flint_malloc(
                                 num_threads * sizeof(_degrees_si_arg_struct));
    start = 0;
    for (i = 0; i < num_threads; i++)
    {
        args[i].degs = degs_array + i*mctx->nvars;
        args[i].start = poly_exps + N*start;
        stop = len*(i+1)/num_threads;
        stop = FLINT_MAX(stop, start);
        stop = FLINT_MIN(stop, len);
        args[i].length = stop - start;
        args[i].bits = bits;
        args[i].mctx = mctx;
        start = stop;
    }

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wake(global_thread_pool,
                                  handles[i], 0, _worker_degrees_si, args + i);
    }

    i = num_handles;
    mpoly_degrees_si(user_degs, args[i].start, args[i].length, bits, mctx);

    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
        for (j = 0; j < mctx->nvars; j++)
        {
            user_degs[j] = FLINT_MAX(user_degs[j], args[i].degs[j]);
        }
    }

    flint_free(degs_array);
    flint_free(args);
}


void mpoly_degrees_ffmpz(
    fmpz * user_degs,
    const ulong * poly_exps,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * max_fields;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            fmpz_set_si(user_degs + i, -WORD(1));
        return;
    }

    TMP_START;

    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields + i);

    mpoly_max_fields_fmpz(max_fields, poly_exps, len, bits, mctx);

    mpoly_get_monomial_ffmpz_unpacked_ffmpz(user_degs, max_fields, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields + i);

    TMP_END;
}

void mpoly_degrees_pfmpz(
    fmpz ** user_degs,
    const ulong * poly_exps,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx)
{
    slong i;
    fmpz * max_fields;
    TMP_INIT;

    if (len == 0)
    {
        for (i = 0; i < mctx->nvars; i++)
            fmpz_set_si(user_degs[i], -WORD(1));
        return;
    }

    TMP_START;

    max_fields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));
    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(max_fields + i);

    mpoly_max_fields_fmpz(max_fields, poly_exps, len, bits, mctx);

    mpoly_get_monomial_pfmpz_unpacked_ffmpz(user_degs, max_fields, mctx);

    for (i = 0; i < mctx->nfields; i++)
        fmpz_clear(max_fields + i);

    TMP_END;
}
