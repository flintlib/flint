/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_urandom(arf_t x, flint_rand_t state, slong bits, arf_rnd_t rnd)
{
    slong prec = bits;
    fmpz_t n;
    fmpz_t t;

    prec += 128;

    fmpz_init(n);
    fmpz_one(n);
    fmpz_mul_2exp(n, n, (ulong) prec);

    fmpz_init(t);
    fmpz_randm(t, state, n);

    arf_set_round_fmpz(x, t, bits, rnd);
    arf_mul_2exp_si(x, x, -prec);

    fmpz_clear(n);
    fmpz_clear(t);
}

