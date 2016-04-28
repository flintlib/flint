/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

slong
fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)
{
    fmpz_t p, q;
    slong i;

    fmpz_init(p);
    fmpz_init(q);

    fmpz_set(p, fmpq_numref(x));
    fmpz_set(q, fmpq_denref(x));

    for (i = 0; i < n && !fmpz_is_zero(q); i++)
    {
        fmpz_fdiv_qr(c + i, p, p, q);
        fmpz_swap(p, q);
    }

    fmpz_set(fmpq_numref(rem), q);
    fmpz_set(fmpq_denref(rem), p);
    fmpq_canonicalise(rem);

    fmpz_clear(p);
    fmpz_clear(q);

    return i;
}
