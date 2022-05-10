/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FQ_NMOD_INLINES_C

#include "fq_nmod.h"

void __fq_nmod_ctx_prime(fmpz_t p, fq_nmod_ctx_t ctx)
{
   fmpz_set(p, fq_nmod_ctx_prime(ctx));
}
