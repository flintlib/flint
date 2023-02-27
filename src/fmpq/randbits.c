/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)
{
    fmpz_randbits(num, state, bits);

    do {
        fmpz_randbits(den, state, bits);
    } while (fmpz_is_zero(den));

    _fmpq_canonicalise(num, den);
}

void fmpq_randbits(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
{
    _fmpq_randbits(fmpq_numref(res), fmpq_denref(res), state, bits);
}
