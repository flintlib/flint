/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            slong r)
{
    fmpz_t u;

    if (r == 0)
    {
        fmpz_zero(rnum);
        fmpz_one(rden);
    } else if (r == 1)
    {
        fmpz_set(rnum, p);
	fmpz_set(rden, q);
    } else if (r == -WORD(1))
    {
        fmpz_neg(rnum, p);
        fmpz_set(rden, q);
    } else
    {
        fmpz_init_set_si(u, r);

	fmpz_gcd(u, u, q);

        if (r == WORD_MIN && !fmpz_fits_si(u)) /* deal with SLONG_MIN */
        {
            fmpz_t fr;
            fmpz_init_set_si(fr, r);
	    fmpz_tdiv_q(fr, fr, u);
	    r = fmpz_get_si(fr);
	    fmpz_tdiv_q(rden, q, u);
	} else if (!fmpz_is_one(u))
        {
            r /= fmpz_get_si(u);
            fmpz_tdiv_q(rden, q, u);
	} else
            fmpz_set(rden, q);

        fmpz_mul_si(rnum, p, r);

	fmpz_clear(u);
    }
}

void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c)
{
    _fmpq_mul_si(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
