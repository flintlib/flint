/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"

void
_fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)
{
    mp_limb_t x = n_randlimb(state);

    fmpz_randtest(num, state, bits);

    if (bits == 1)
    {
        fmpz_one(den);
        return;
    }

    fmpz_randtest_not_zero(den, state, bits);

    switch (x % 16)
    {
        case 0:
            fmpz_set_si(num, 1);
            break;
        case 1:
            fmpz_set_si(num, -1);
            break;
        case 2:
            fmpz_set_si(num, 2);
            break;
        case 3:
            fmpz_set_si(num, -2);
            break;
    }

    switch ((x / 16) % 16)
    {
        case 0:
            fmpz_set_si(den, 1);
            break;
        case 2:
            fmpz_set_si(den, 2);
            break;
    }

    _fmpq_canonicalise(num, den);
}

void fmpq_randtest(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
{
    _fmpq_randtest(fmpq_numref(res), fmpq_denref(res), state, bits);
}

void fmpq_randtest_not_zero(fmpq_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    if (bits == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_randtest_not_zero). bits == 0.\n");
    }

    do {
        fmpq_randtest(f, state, bits);
    } while (fmpq_is_zero(f));
}
