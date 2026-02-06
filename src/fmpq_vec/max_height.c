/*
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"

void
_fmpq_vec_max_height(fmpz_t height, const fmpq * vec, slong len)
{
    fmpz_t h;
    slong i;

    fmpz_init(h);
    fmpz_zero(height);

    for (i = 0; i < len; i++)
    {
        fmpq_height(h, vec + i);
        if (fmpz_cmpabs(h, height) > 0)
            fmpz_swap(height, h);
    }

    fmpz_clear(h);
}
