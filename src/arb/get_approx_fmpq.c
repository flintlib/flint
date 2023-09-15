/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb.h"

int arb_get_approx_fmpq(fmpq_t y, const arb_t x, slong prec)
{
    int res;
    mpz_t n, d;
    fmpz_t c;
    arb_t z;

    mpz_init(n);
    mpz_init(d);
    fmpz_init(c);
    arb_init(z);

    res = arf_get_approx_fmpq(y, arb_midref(x), prec);

    if (res)
    {
        fmpq_get_mpz_frac(n, d, y);
        fmpz_set_mpz(c, d);
        arb_mul_fmpz(z, x, c, prec);
        res = (mag_cmp_2exp_si(arb_radref(z), - prec/8) < 0);
    }

    mpz_clear(n);
    mpz_clear(d);
    fmpz_clear(c);
    arb_clear(z);
    return res;
}
