/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
#include <math.h>
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "long_extras.h"

typedef struct
{
    fmpz * res;
    const fmpz * vec;
    const fmpz_mod_ctx_struct * ctx;
}
work_t;

static void
worker(slong i, void * args)
{
    work_t * w = (work_t *) args;
    fmpz_mod_set_fmpz(w->res + i, w->vec + i, w->ctx);
}

static void
_fmpz_mod_vec_set_fmpz_vec_threaded(fmpz * A, const fmpz * B, slong len,
                                                     const fmpz_mod_ctx_t ctx)
{
    work_t work[1];

    work->res = A;
    work->vec = B;
    work->ctx = ctx;

    flint_parallel_do(worker, work, len, 0, FLINT_PARALLEL_UNIFORM);
}

void
_fmpz_mod_vec_set_fmpz_vec(fmpz * A, const fmpz * B, slong len,
                                                     const fmpz_mod_ctx_t ctx)
{
    if (len >= 2)
    {
        slong bits = fmpz_bits(fmpz_mod_ctx_modulus(ctx));

        if ((len >= 10000 ||
            (bits >= 20000 && fabs((double) _fmpz_vec_max_bits(B, len)) >= 20000) ||
                (len * (double) bits >= 400000 && len * fabs((double)_fmpz_vec_max_bits(B, len)) >= 400000)) &&
            flint_get_num_threads() >= 2)
        {
            _fmpz_mod_vec_set_fmpz_vec_threaded(A, B, len, ctx);
            return;
        }
    }

    for (len--; len >= 0; len--)
        fmpz_mod_set_fmpz(A + len, B + len, ctx);
}
