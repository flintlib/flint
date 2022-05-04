/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fq_nmod.h"

void
fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx)
{
    fmpz_set(f, fq_nmod_ctx_prime(ctx));
    fmpz_pow_ui(f, f, fq_nmod_ctx_degree(ctx));
}
