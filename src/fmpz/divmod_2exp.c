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

/* Hensel division: sets rop to the unique q in [0, 2^N) with
   q b == a (mod 2^N). Requires b odd. */
void
fmpz_divmod_2exp(fmpz_t rop, const fmpz_t a, const fmpz_t b, flint_bitcnt_t N)
{
    fmpz_t ta, tb;
    mp_size_t n, an, bn;
    mpz_ptr ma, mb, mr;
    mp_ptr rd;

    if (!fmpz_is_odd(b))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_divmod_2exp). b is even.\n");
    }

    if (N == 0)
    {
        fmpz_zero(rop);
        return;
    }

    n = (N + FLINT_BITS - 1) / FLINT_BITS;

    fmpz_init(ta);
    fmpz_init(tb);
    fmpz_fdiv_r_2exp(ta, a, FLINT_BITS * n);
    fmpz_fdiv_r_2exp(tb, b, FLINT_BITS * n);

    if (fmpz_is_zero(ta))
    {
        fmpz_zero(rop);
    }
    else
    {
        ma = _fmpz_promote_val(ta);
        mb = _fmpz_promote_val(tb);
        an = ma->_mp_size;
        bn = mb->_mp_size;

        mr = _fmpz_promote(rop);
        rd = FLINT_MPZ_REALLOC(mr, n);
        flint_mpn_bdiv_q(rd, ma->_mp_d, an, mb->_mp_d, bn, n);
        while (n > 0 && rd[n - 1] == 0)
            n--;
        mr->_mp_size = n;
        _fmpz_demote_val(rop);
        fmpz_fdiv_r_2exp(rop, rop, N);
    }

    fmpz_clear(ta);
    fmpz_clear(tb);
}
