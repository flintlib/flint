/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "bernoulli.h"

void
_bernoulli_fmpq_ui(fmpz_t num, fmpz_t den, ulong n)
{
    if (n < (ulong) bernoulli_cache_num)
    {
        fmpz_set(num, fmpq_numref(bernoulli_cache + n));
        fmpz_set(den, fmpq_denref(bernoulli_cache + n));
    }
    else if (n < 18000 || n % 2 == 1)
    {
        _bernoulli_fmpq_ui_zeta(num, den, n);
    }
    else
    {
        _bernoulli_fmpq_ui_multi_mod(num, den, n, -1.0);
    }
}

void
bernoulli_fmpq_ui(fmpq_t b, ulong n)
{
    _bernoulli_fmpq_ui(fmpq_numref(b), fmpq_denref(b), n);
}

