/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Returns 1 if gcd(x * y, n) == 1; otherwise returns 0.
*/
int
aprcl_is_mul_coprime_ui_ui(ulong x, ulong y, const fmpz_t n)
{
    ulong a, rem;
    int result = 0;

    rem = fmpz_tdiv_ui(n, x);
    a = n_gcd(x, rem);              /* a = gcd(x, n % x) */
    if (a == 1)
    {
       rem = fmpz_tdiv_ui(n, y);
       result = n_gcd(y, rem) == 1;     /* result = gcd(y, m % y) */
    }

    return result;
}


/*
    Returns 1 if gcd(x * y, n) == 1; otherwise returns 0.
*/
int
aprcl_is_mul_coprime_ui_fmpz(ulong x, const fmpz_t y, const fmpz_t n)
{
    int is_coprime = 0;
    ulong a, rem;
    fmpz_t result;

    fmpz_init(result);

    rem = fmpz_tdiv_ui(n, x);
    a = n_gcd(x, rem);              /* a = gcd(x, n % x) */
    if (a == 1)
    {
       fmpz_fdiv_r(result, n, y);      /* result = n % y */
       fmpz_gcd(result, result, y);    /* result = gcd(y, n % y) */

       is_coprime = fmpz_is_one(result);
    }

    fmpz_clear(result);

    return is_coprime;
}
