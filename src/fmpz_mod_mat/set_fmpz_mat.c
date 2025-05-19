/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"

typedef struct {
    fmpz * M;
    slong Mstride;
    const fmpz * A;
    slong Astride;
    slong c;
    const fmpz_mod_ctx_struct * ctx;
} _worker_arg;

static void _red_worker(slong i, void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    _fmpz_mod_vec_set_fmpz_vec(arg->M + i * arg->Mstride, arg->A + i * arg->Astride, arg->c, arg->ctx);
}

static void
_fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t M, const fmpz_mat_t A, const fmpz_mod_ctx_t ctx)
{
    slong i, r, c;

    r = fmpz_mod_mat_nrows(M, ctx);
    c = fmpz_mod_mat_ncols(M, ctx);

    for (i = 0; i < r; i++)
        _fmpz_mod_vec_set_fmpz_vec(fmpz_mod_mat_row(M, i), fmpz_mat_row(A, i), c, ctx);
}

void fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t M, const fmpz_mat_t A, const fmpz_mod_ctx_t ctx)
{
    slong r, c;
    slong limit;

    r = fmpz_mod_mat_nrows(M, ctx);
    c = fmpz_mod_mat_ncols(M, ctx);

    /* limit on threads */
    limit = fmpz_size(ctx->n) + r + c;
    limit = limit < 64 ? 0 : (limit - 64)/64;
    limit = FLINT_MIN(limit, r);

    if (limit < 2)
    {
        _fmpz_mod_mat_set_fmpz_mat(M, A, ctx);
    }
    else
    {
        _worker_arg arg;

        arg.M = M->entries;
        arg.Mstride = M->stride;
        arg.A = A->entries;
        arg.Astride = A->stride;
        arg.c = c;
        arg.ctx = ctx;

        flint_parallel_do(_red_worker, &arg, r, limit, FLINT_PARALLEL_UNIFORM);
    }
}
