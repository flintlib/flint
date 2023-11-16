/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

/* ANDs used without zeroing */
void
_fmpz_mpoly_used_vars(int * used, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)
{
    slong i, n;
    fmpz * max;
    TMP_INIT;

    n = ctx->minfo->nvars;

    if (fmpz_mpoly_is_fmpz(f, ctx))
        return;

    /* The polynomial is not constant and has only one variable. */
    if (n == 1)
    {
        used[0] = 1;
        return;
    }

    TMP_START;
    max = TMP_ALLOC(sizeof(fmpz) * n);
    for (i = 0; i < n; i++)
        fmpz_init(max + i);

    mpoly_degrees_ffmpz(max, f->exps, f->length, f->bits, ctx->minfo);

    for (i = 0; i < n; i++)
        used[i] |= !fmpz_is_zero(max + i);

    for (i = 0; i < n; i++)
        fmpz_clear(max + i);

    TMP_END;
}

void
fmpz_mpoly_q_used_vars(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
{
    slong i, n;

    n = ctx->minfo->nvars;

    for (i = 0; i < n; i++)
        used[i] = 0;

    _fmpz_mpoly_used_vars(used, fmpz_mpoly_q_numref(f), ctx);
    _fmpz_mpoly_used_vars(used, fmpz_mpoly_q_denref(f), ctx);
}

void
fmpz_mpoly_q_used_vars_num(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
{
    slong i, n;

    n = ctx->minfo->nvars;

    for (i = 0; i < n; i++)
        used[i] = 0;

    _fmpz_mpoly_used_vars(used, fmpz_mpoly_q_numref(f), ctx);
}

void
fmpz_mpoly_q_used_vars_den(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
{
    slong i, n;

    n = ctx->minfo->nvars;

    for (i = 0; i < n; i++)
        used[i] = 0;

    _fmpz_mpoly_used_vars(used, fmpz_mpoly_q_denref(f), ctx);
}
