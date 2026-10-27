/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_vec.h"

void
_fmpq_vec_get_fmpz_vec_fmpz(fmpz* num, fmpz_t den, const fmpq * a, slong len)
{
    slong i;

    if (len < 1)
    {
        fmpz_one(den);
        return;
    }
    else if (len == 1)
    {
        fmpz_set(num + 0, fmpq_numref(a + 0));
        fmpz_set(den, fmpq_denref(a + 0));
        return;
    }

    /* the denominators are every second fmpz */
    _fmpz_vec_lcm_strided(den, fmpq_denref(a), 2, len, NULL);

    if (fmpz_is_one(den))
    {
        for (i = 0; i < len; i++)
            fmpz_set(num + i, fmpq_numref(a + i));
    }
    else
    {
        _fmpq_vec_get_fmpz_vec_scaled(num, 1, a, 1, len, den);
    }
}
