/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            slong r)
{
    fmpz_t u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && r >= COEFF_MIN && r <= COEFF_MAX)
    {
        _fmpq_add_small(rnum, rden, *p, *q, r, 1);
        return;
    }

    /* both are integers */
    if (fmpz_is_one(q))
    {
        if (r >= 0)
           fmpz_add_ui(rnum, p, r);
        else
           fmpz_sub_ui(rnum, p, -r);

        fmpz_set(rden, q);
        
        return;
    }

    /*
    We want to compute p/q + r/1 where the inputs are already
    in canonical form.

    Note (p + q*r, q) is in canonical form.

    */
    
    fmpz_init(u);

    fmpz_mul_si(u, q, r);
    fmpz_add(rnum, p, u);
    fmpz_set(rden, q);

    fmpz_clear(u);
}

void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c)
{
    _fmpq_add_si(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
