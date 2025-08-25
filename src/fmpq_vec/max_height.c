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
    slong i;
    fmpz * Lheights;
    Lheights = _fmpz_vec_init(len);

    for (i = 0; i < len; i++)
        fmpq_height(Lheights + i, vec + i);

    _fmpz_vec_height(height, Lheights, len);
}
