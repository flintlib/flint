/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Memory management. */

void
mp_real_init(mp_real_t x)
{
    x->d = NULL;
    x->alloc = 0;
    x->size = 0;
    x->negative = 0;
    x->exp = 0;
    x->err = 0;
}

void
mp_real_clear(mp_real_t x)
{
    flint_free(x->d);
}

void
_mp_real_grow(mp_real_t x, slong k)
{
    slong newalloc = FLINT_MAX(k, x->alloc + x->alloc / 2);
    x->d = flint_realloc(x->d, newalloc * sizeof(ulong));
    x->alloc = newalloc;
}

void
mp_real_zero(mp_real_t x)
{
    x->size = 0;
    x->negative = 0;
    x->exp = 0;
    x->err = 0;
}

void
mp_real_swap(mp_real_t x, mp_real_t y)
{
    FLINT_SWAP(mp_real_struct, *x, *y);
}
