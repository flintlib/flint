/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"

/*
    If op is a square modulo 2^N (for N >= 3, iff op == 1 mod 8; for N = 2
    iff op == 1 mod 4; for N = 1 iff op is odd), sets rop to a square root
    modulo 2^N (the one congruent to 1 modulo 4, reduced modulo 2^N) and
    returns 1; otherwise returns 0. Even op is not supported (returns 0).
*/
int
fmpz_sqrtmod_2exp(fmpz_t rop, const fmpz_t op, flint_bitcnt_t N)
{
    fmpz_t t;
    mp_size_t n, xn;
    mpz_ptr mt, mr;
    mp_ptr rd;
    ulong low = fmpz_fdiv_ui(op, 8);

    if (N == 0)
    {
        fmpz_zero(rop);
        return 1;
    }

    if ((low & 1) == 0)
        return 0;

    if (N <= 3)
    {
        if ((N == 2 && (low & 3) != 1) || (N == 3 && low != 1))
            return 0;
        fmpz_one(rop);
        return 1;
    }

    if (low != 1)
        return 0;

    n = (N + FLINT_BITS - 1) / FLINT_BITS;

    fmpz_init(t);
    fmpz_fdiv_r_2exp(t, op, FLINT_BITS * n);
    mt = _fmpz_promote_val(t);
    xn = mt->_mp_size;

    mr = _fmpz_promote(rop);
    rd = FLINT_MPZ_REALLOC(mr, n);
    flint_mpn_bsqrt(rd, mt->_mp_d, xn, n);
    while (n > 0 && rd[n - 1] == 0)
        n--;
    mr->_mp_size = n;
    _fmpz_demote_val(rop);

    fmpz_fdiv_r_2exp(rop, rop, N);
    fmpz_clear(t);
    return 1;
}
