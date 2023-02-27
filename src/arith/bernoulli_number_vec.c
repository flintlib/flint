/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void _arith_bernoulli_number_vec(fmpz * num, fmpz * den, slong n)
{
    if (n < 700)
        _arith_bernoulli_number_vec_recursive(num, den, n);
    else if (n < 3900)
        _arith_bernoulli_number_vec_zeta(num, den, n);
    else
        _arith_bernoulli_number_vec_multi_mod(num, den, n);
}

void arith_bernoulli_number_vec(fmpq * x, slong n)
{
    fmpz * num, * den;
    slong i;

    if (n <= 0)
        return;

    num = _fmpz_vec_init(n * 2);
    den = num + n;

    _arith_bernoulli_number_vec(num, den, n);

    for (i = 0; i < n; i++)
    {
        fmpz_swap(num + i, fmpq_numref(x + i));
        fmpz_swap(den + i, fmpq_denref(x + i));
    }

    _fmpz_vec_clear(num, n * 2);
}

