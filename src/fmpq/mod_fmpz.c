/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int
_fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod)
{
    int result;
    fmpz_t tmp;

    fmpz_init(tmp);
    result = fmpz_invmod(tmp, den, mod);
    fmpz_mul(tmp, tmp, num);
    fmpz_mod(res, tmp, mod);
    fmpz_clear(tmp);

    return result;
}

int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)
{
    return _fmpq_mod_fmpz(res, fmpq_numref(x), fmpq_denref(x), mod);
}

