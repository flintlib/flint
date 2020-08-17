/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

/* TODO: optimize */
void
fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx)
{
    fmpz_t exp;
    fmpz_init_set_ui(exp, e);
    fq_pow(rop, op, exp, ctx);
    fmpz_clear(exp);
}
