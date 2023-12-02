/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden,
                  const fmpz_t opnum, const fmpz_t opden, slong e)
{
    if (e >= 0)
    {
        fmpz_pow_ui(rnum, opnum, e);
        fmpz_pow_ui(rden, opden, e);
    }
    else
    {
        if (rnum == opnum)
        {
            fmpz t;

            fmpz_pow_ui(rnum, opnum, -e);
            fmpz_pow_ui(rden, opden, -e);

            t     = *rnum;
            *rnum = *rden;
            *rden = t;
        }
        else
        {
            fmpz_pow_ui(rden, opnum, -e);
            fmpz_pow_ui(rnum, opden, -e);
        }

        if (fmpz_sgn(rden) < 0)
        {
            fmpz_neg(rnum, rnum);
            fmpz_neg(rden, rden);
        }
    }
}

void fmpq_pow_si(fmpq_t rop, const fmpq_t op, slong e)
{
    if (e < 0 && fmpz_is_zero(fmpq_numref(op)))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_pow_si). Division by zero.\n");
    }

    _fmpq_pow_si(fmpq_numref(rop), fmpq_denref(rop),
                 fmpq_numref(op), fmpq_denref(op), e);
}

