/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx,
                                             slong nvars, const ordering_t ord)
{
   ctx->n = (ord == ORD_DEGLEX || ord == ORD_DEGREVLEX) ? nvars + 1 : nvars;
   ctx->ord = ord;
}
