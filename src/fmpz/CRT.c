/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
_fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2,
                   const fmpz_t m2, const fmpz_t m1m2, fmpz_t c, int sign)
{
    fmpz_t r1normal, tmp, r1mod, s;

    fmpz_init(tmp);
    fmpz_init(r1mod);
    fmpz_init(s);

    /* FIXME: assume r1 moved to [0, m1); add tests for this */
    if (fmpz_sgn(r1) < 0)
    {
        fmpz_init(r1normal);
        fmpz_add(r1normal, r1, m1);
    }
    else
    {
        *r1normal = *r1;
    }

    fmpz_mod(r1mod, r1normal, m2);
    fmpz_sub(s, r2, r1mod);
    if (fmpz_sgn(s) < 0)
       fmpz_add(s, s, m2);
    fmpz_mul(s, s, c);
    fmpz_mod(s, s, m2);
    fmpz_mul(tmp, m1, s);
    fmpz_add(tmp, tmp, r1normal);

    if (fmpz_sgn(r1) < 0)
        fmpz_clear(r1normal);

    if (sign)
    {
        fmpz_sub(out, tmp, m1m2);
        if (fmpz_cmpabs(tmp, out) <= 0)
            fmpz_set(out, tmp);
    }
    else
    {
        fmpz_set(out, tmp);
    }

    fmpz_clear(tmp);
    fmpz_clear(r1mod);
    fmpz_clear(s);
}

void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
    const fmpz_t r2, const fmpz_t m2, int sign)
{
    fmpz_t m1m2, c;

    fmpz_init(c);

    fmpz_mod(c, m1, m2);
    if (!fmpz_invmod(c, c, m2))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_CRT). m1 not invertible modulo m2.\n");
    }

    fmpz_init(m1m2);
    fmpz_mul(m1m2, m1, m2);

    _fmpz_CRT(out, r1, m1, r2, m2, m1m2, c, sign);

    fmpz_clear(m1m2);
    fmpz_clear(c);
}
