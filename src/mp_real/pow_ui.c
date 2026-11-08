/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mp_real.h"

/* left-to-right binary powering, each product at n limbs */
void
mp_real_pow_ui(mp_real_t res, const mp_real_t x, ulong e, slong n)
{
    mp_real_t t;
    const mp_real_struct * b = x;
    slong i;

    if (e == 0)
    {
        mp_real_set_ui(res, 1);
        return;
    }

    mp_real_init(t);
    if (res == x)
    {
        mp_real_set(t, x);
        b = t;
    }

    mp_real_set(res, b);
    for (i = FLINT_BIT_COUNT(e) - 2; i >= 0; i--)
    {
        mp_real_mul(res, res, res, n);
        if ((e >> i) & 1)
            mp_real_mul(res, res, b, n);
    }

    mp_real_clear(t);
}
