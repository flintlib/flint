/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


void _fmpq_set_cfrac_divconquer(_fmpz_mat22_t M, const fmpz * c, slong n)
{
    _fmpz_mat22_one(M);
    if (n < 32)
    {
        slong i;
        for (i = 0; i < n; i++)
            _fmpz_mat22_rmul_elem(M, c + i);
    }
    else
    {
        slong m = n / 2;
        _fmpz_mat22_t N;
        _fmpz_mat22_init(N);
        _fmpq_set_cfrac_divconquer(M, c, m);
        _fmpq_set_cfrac_divconquer(N, c + m, n - m);
        _fmpz_mat22_rmul(M, N);
        _fmpz_mat22_clear(N);
    }
}


void fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n)
{
    _fmpz_mat22_t M;
    _fmpz_mat22_init(M);
    _fmpq_set_cfrac_divconquer(M, c, n);
    fmpz_swap(fmpq_numref(x), M->_11);
    fmpz_swap(fmpq_denref(x), M->_21);
    _fmpz_mat22_clear(M);
    FLINT_ASSERT(n <= 0 || fmpq_is_canonical(x));
}

