/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

/* get the (unpacked) fields for the generator of index var */
void mpoly_gen_fields_ui(ulong * gexp, slong var, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong nfields = mctx->nfields;
    int deg = mctx->deg;
    int rev = mctx->rev;
    slong i;

    for (i = 0; i < nfields; i++)
        gexp[i] = WORD(0);

    gexp[rev ? var : nvars - 1 - var] = WORD(1);
    if (deg)
        gexp[nvars] = WORD(1);

}

void mpoly_gen_fields_fmpz(fmpz * gexp, slong var, const mpoly_ctx_t mctx)
{
    slong nvars = mctx->nvars;
    slong nfields = mctx->nfields;
    int deg = mctx->deg;
    int rev = mctx->rev;
    slong i;

    for (i = 0; i < nfields; i++)
        fmpz_zero(gexp + i);

    fmpz_one(gexp + (rev ? var : nvars - 1 - var));
    if (deg)
        fmpz_one(gexp + nvars);
}

