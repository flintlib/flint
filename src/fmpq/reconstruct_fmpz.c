/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"


int _fmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m)
{
    fmpz_t N;
    int result;

    fmpz_init(N);
    fmpz_fdiv_q_2exp(N, m, 1);
    if (fmpz_is_even(m))
        fmpz_sub_ui(N, N, 1);
    fmpz_sqrt(N, N);
    result = _fmpq_reconstruct_fmpz_2(n, d, a, m, N, N);
    fmpz_clear(N);

    return result;
}

int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
{
    return _fmpq_reconstruct_fmpz(fmpq_numref(res), fmpq_denref(res), a, m);
}

