/*
    Copyright (C) 2011, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_prod(fmpz_t res, const fmpz * vec, slong len)
{
    if (len <= 1)
    {
        if (len == 1)
            fmpz_set(res, vec);
        else
            fmpz_one(res);
    }
    else if (len <= 3)
    {
        slong i;

        fmpz_mul(res, vec, vec + 1);
        for (i = 2; i < len; i++)
            fmpz_mul(res, res, vec + i);
    }
    else
    {
        slong m = len / 2;
        fmpz_t tmp;
        fmpz_init(tmp);
        _fmpz_vec_prod(res, vec, m);
        _fmpz_vec_prod(tmp, vec + m, len - m);
        fmpz_mul(res, res, tmp);
        fmpz_clear(tmp);
    }
}

void
_fmpz_ui_vec_prod(fmpz_t res, nn_srcptr vec, slong len)
{
    if (!len)
        fmpz_one(res);
    else if (len < 16)
    {
        slong i;

        fmpz_set_ui(res, vec[0]);
        for (i = 1; i < len; i++)
            fmpz_mul_ui(res, res, vec[i]);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_ui_vec_prod(res, vec, len / 2);
        _fmpz_ui_vec_prod(t, vec + len / 2, len - len / 2);
        fmpz_mul(res, res, t);
        fmpz_clear(t);
    }
}
