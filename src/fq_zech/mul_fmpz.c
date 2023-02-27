/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void
fq_zech_mul_fmpz(fq_zech_t rop, const fq_zech_t op, const fmpz_t x,
                 const fq_zech_ctx_t ctx)
{
    mp_limb_t ux;
    fmpz_t y;

    fmpz_init(y);
    fmpz_mod_ui(y, x, ctx->p);

    ux = fmpz_get_ui(y);

    fmpz_clear(y);

    fq_zech_mul_ui(rop, op, ux, ctx);
}
