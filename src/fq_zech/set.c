/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fq_zech.h"

void
fq_zech_set_si(fq_zech_t rop, const slong x, const fq_zech_ctx_t ctx)
{
    fmpz_t xx;
    fmpz_init_set_si(xx, x);
    fq_zech_set_fmpz(rop, xx, ctx);
    fmpz_clear(xx);
}

void
fq_zech_set_ui(fq_zech_t rop, const ulong x, const fq_zech_ctx_t ctx)
{
    fmpz_t xx;
    fmpz_init_set_ui(xx, x);
    fq_zech_set_fmpz(rop, xx, ctx);
    fmpz_clear(xx);
}
