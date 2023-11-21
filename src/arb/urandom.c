/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_urandom(arb_t x, flint_rand_t state, slong bits)
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

    arb_set_round_fmpz(x, t, bits);
    arb_mul_2exp_si(x, x, -prec);

    fmpz_clear(n);
    fmpz_clear(t);
}

