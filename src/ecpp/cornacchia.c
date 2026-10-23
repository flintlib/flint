/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "ecpp.h"

/*
    Modified Cornacchia (Cohen, Algorithm 1.5.3): given a discriminant
    D < 0, an odd n with gcd(n, D) = 1 and a square root sqrtD of D
    modulo n, finds t, v >= 0 with t^2 + |D| v^2 = 4n. Returns 1 on
    success, 0 if there is no solution (or n is composite and the
    algorithm fails).
*/
int
ecpp_cornacchia(fmpz_t t, fmpz_t v, const fmpz_t n, slong D, const fmpz_t sqrtD)
{
    fmpz_t x0, a, b, l, r, absD;
    int result = 0;

    fmpz_init(x0);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(l);
    fmpz_init(r);
    fmpz_init_set_si(absD, -D);

    /* x0 = sqrt(D) mod n, with x0 = D mod 2, 0 < x0 < 2n */
    fmpz_mod(x0, sqrtD, n);
    if (fmpz_is_odd(x0) != (FLINT_ABS(D) % 2 == 1))
        fmpz_sub(x0, n, x0);
    /* now x0 = D mod 2 provided n odd; x0 is a root of x^2 = D mod 4n */

    /* Euclid on (2n, x0) until the remainder drops below 2 sqrt(n) */
    fmpz_mul_2exp(a, n, 1);
    fmpz_set(b, x0);
    fmpz_mul_2exp(l, n, 2);
    fmpz_sqrt(l, l);            /* l = floor(sqrt(4n)) */

    /* Lehmer partial gcd: (a, b) become the last two remainders, b <= l */
    if (fmpz_cmp(b, l) > 0)
    {
        fmpz_t co1, co2;
        fmpz_init(co1);
        fmpz_init(co2);
        fmpz_xgcd_partial(co2, co1, a, b, l);
        fmpz_clear(co1);
        fmpz_clear(co2);
    }

    /* b^2 + |D| v^2 = 4n ? */
    fmpz_mul(r, b, b);
    fmpz_mul_2exp(l, n, 2);
    fmpz_sub(r, l, r);          /* r = 4n - b^2 */
    if (fmpz_sgn(r) > 0 && fmpz_divisible(r, absD))
    {
        fmpz_divexact(r, r, absD);
        if (fmpz_is_square(r))
        {
            fmpz_sqrt(v, r);
            fmpz_set(t, b);
            result = 1;
        }
    }

    fmpz_clear(x0);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(l);
    fmpz_clear(r);
    fmpz_clear(absD);

    return result;
}
