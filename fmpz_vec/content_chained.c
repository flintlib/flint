/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"

void
_fmpz_vec_content_chained(fmpz_t res, const fmpz * vec, slong len, const fmpz_t in)
{
    while (len > 0 && fmpz_is_zero(vec + 0))
    {
        len--;
        vec++;
    }

    while (len > 1 && fmpz_is_zero(vec + len - 1))
        len--;

    if (len == 0)
    {
        fmpz_abs(res, in);
        return;
    }

    if (len == 1)
    {
        fmpz_gcd(res, vec + 0, in);
        return;
    }

    if (fmpz_is_pm1(in) || fmpz_is_pm1(vec + 0) || fmpz_is_pm1(vec + len - 1))
    {
        fmpz_one(res);
        return;
    }

    fmpz_gcd3(res, vec + 0, vec + len - 1, in);
    vec++;
    len -= 2;

    while (len >= 2 && !fmpz_is_one(res))
    {
        fmpz_gcd3(res, vec + 0, vec + len - 1, res);
        vec++;
        len -= 2;
    }

    if (len != 0 && !fmpz_is_one(res))
        fmpz_gcd(res, res, vec + 0);
}
