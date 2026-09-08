/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

/*
    Checked exact square root: if a is a perfect square, sets s to its
    nonnegative square root and returns 1; otherwise returns 0 with s
    undefined. Nonsquares are usually rejected cheaply by quadratic residue
    tests inside flint_mpn_sqrt.
*/
int
fmpz_perfect_sqrt(fmpz_t s, const fmpz_t a)
{
    fmpz c = *a;

    if (!COEFF_IS_MPZ(c))
    {
        ulong r, t;

        if (c < 0)
            return 0;

        t = n_sqrtrem(&r, (ulong) c);
        if (r != 0)
            return 0;
        fmpz_set_ui(s, t);
        return 1;
    }
    else
    {
        mpz_srcptr ma = COEFF_TO_PTR(c);
        mp_size_t an = ma->_mp_size, sn;
        mpz_ptr ms;
        mp_ptr sd;
        int square;

        if (an < 0)
            return 0;

        /* cheap rejection before touching s: squares modulo 256 */
        {
            static const unsigned char sq256[32] = {0x13, 0x02, 0x03, 0x02, 0x12, 0x02, 0x02, 0x02, 0x13, 0x02, 0x02, 0x02, 0x12, 0x02, 0x02, 0x02, 0x12, 0x02, 0x03, 0x02, 0x12, 0x02, 0x02, 0x02, 0x12, 0x02, 0x02, 0x02, 0x12, 0x02, 0x02, 0x02};
            mp_limb_t x = ma->_mp_d[0] & 255;
            if (!((sq256[x >> 3] >> (x & 7)) & 1))
                return 0;
        }

        if (s == a)
        {
            fmpz_t t;
            fmpz_init(t);
            square = fmpz_perfect_sqrt(t, a);
            fmpz_swap(s, t);
            fmpz_clear(t);
            return square;
        }

        sn = (an + 1) / 2;
        ms = _fmpz_promote(s);
        sd = FLINT_MPZ_REALLOC(ms, sn);
        square = flint_mpn_sqrt(sd, ma->_mp_d, an);
        ms->_mp_size = square ? sn : 0;
        _fmpz_demote_val(s);
        return square;
    }
}
