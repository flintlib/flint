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

/* rop = op^(-1) mod 2^N for odd op (via flint_mpn_binv) */
void
fmpz_invmod_2exp(fmpz_t rop, const fmpz_t op, flint_bitcnt_t N)
{
    fmpz_t t;
    mp_size_t n, xn;
    mpz_ptr mt, mr;
    mp_ptr rd;

    if (!fmpz_is_odd(op))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_invmod_2exp). op is even.\n");
    }

    if (N == 0)
    {
        fmpz_zero(rop);
        return;
    }

    n = (N + FLINT_BITS - 1) / FLINT_BITS;

    if (n == 1)
    {
        ulong x0, r;

        /* op mod 2^FLINT_BITS */
        if (!COEFF_IS_MPZ(*op))
            x0 = (ulong) *op;
        else if (COEFF_TO_PTR(*op)->_mp_size > 0)
            x0 = COEFF_TO_PTR(*op)->_mp_d[0];
        else
            x0 = -COEFF_TO_PTR(*op)->_mp_d[0];

        r = n_binvert(x0);
        if (N < FLINT_BITS)
            r &= (UWORD(1) << N) - 1;
        fmpz_set_ui(rop, r);
        return;
    }

    /* t = op mod 2^(FLINT_BITS n) as a nonnegative multi-limb integer */
    fmpz_init(t);
    fmpz_fdiv_r_2exp(t, op, FLINT_BITS * n);
    mt = _fmpz_promote_val(t);
    xn = mt->_mp_size;

    mr = _fmpz_promote(rop);
    rd = FLINT_MPZ_REALLOC(mr, n);
    flint_mpn_binv(rd, mt->_mp_d, xn, n);
    while (n > 0 && rd[n - 1] == 0)
        n--;
    mr->_mp_size = n;
    _fmpz_demote_val(rop);

    fmpz_fdiv_r_2exp(rop, rop, N);
    fmpz_clear(t);
}
