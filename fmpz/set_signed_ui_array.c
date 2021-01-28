/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

/*
    Given an array of limbs "in" representing a integer mod 2^(FLINT_BITS*n),
    set "out" to the symmetric remainder with the halfway point
    2^(FLINT_BITS*n/2) mapping to -2^(FLINT_BITS*n/2)
*/

void fmpz_set_signed_ui_array(fmpz_t f, const ulong * c_in, slong n)
{
    slong i;
    ulong * c;
    int neg;

    TMP_INIT;

    TMP_START;

    c = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    for (i = 0; i < n; i++)
       c[i] = c_in[i];

    neg = 0 > (slong) c[n - 1];

    if (neg)
       mpn_neg_n(c, c, n);

    while (n > 0 && c[n - 1] == 0)
       n--;

    if (n <= 1)
    {
       fmpz_set_ui(f, c[0]);
    }
    else
    {
       __mpz_struct * mpz = _fmpz_promote(f);

       mpz_realloc2(mpz, n*FLINT_BITS);

       mpn_copyi(mpz->_mp_d, c, n);
       mpz->_mp_size = n;
    }

    if (neg)
       fmpz_neg(f, f);

    TMP_END;
}

