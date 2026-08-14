/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_nmod.h"

void
padic_nmod_ctx_clear(gr_ctx_t ctx)
{
    flint_free(PADIC_NMOD_CTX_POW(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}
