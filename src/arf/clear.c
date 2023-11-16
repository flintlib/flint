/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_clear(arf_t x)
{
    ARF_DEMOTE(x);

    /* fmpz_clear(), but avoids a redundant store */
    if (COEFF_IS_MPZ(ARF_EXP(x)))
        _fmpz_clear_mpz(ARF_EXP(x));
}

arf_ptr
_arf_vec_init(slong n)
{
    slong i;
    arf_ptr v = (arf_ptr) flint_malloc(sizeof(arf_struct) * n);

    for (i = 0; i < n; i++)
        arf_init(v + i);

    return v;
}

void
_arf_vec_clear(arf_ptr x, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arf_clear(x + i);
    flint_free(x);
}
