/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FQ_INLINES_C

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

#include "fq.h"

#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

#include "fmpz.h"

/* TODO: Remove this */
void __fq_ctx_prime(fmpz_t p, fq_ctx_t ctx)
{
    fmpz_set(p, fq_ctx_prime(ctx));
}
