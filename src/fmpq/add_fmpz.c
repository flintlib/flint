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
_fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            const fmpz_t r)
{
    fmpz_t u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && !COEFF_IS_MPZ(*r))
    {
        _fmpq_add_small(rnum, rden, *p, *q, *r, 1);
        return;
    }

    /* both are integers */
    if (fmpz_is_one(q))
    {
        fmpz_add(rnum, p, r);
        
        fmpz_set(rden, q);
        
        return;
    }

    /*
    We want to compute p/q + r/1 where the inputs are already
    in canonical form.

    Note (p + q*r, q) is in canonical form.

    */
    
    fmpz_init(u);

    fmpz_mul(u, q, r);
    fmpz_add(rnum, p, u);
    fmpz_set(rden, q);

    fmpz_clear(u);
}

void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
{
    _fmpq_add_fmpz(fmpq_numref(res), fmpq_denref(res),
              fmpq_numref(op1), fmpq_denref(op1), c);
}
