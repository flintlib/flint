/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq.h"

void
_fmpq_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden, 
                        const fmpz * poly, const fmpz_t den, slong len, 
                        const fmpz_t anum, const fmpz_t aden)
{
    fmpz_t d;
    
    _fmpz_poly_evaluate_fmpq(rnum, rden, poly, len, anum, aden);
    fmpz_mul(rden, rden, den);
    
    fmpz_init(d);
    fmpz_gcd(d, rnum, rden);
    if (!fmpz_is_one(d))
    {
        fmpz_divexact(rnum, rnum, d);
        fmpz_divexact(rden, rden, d);
    }
    fmpz_clear(d);
}

void 
fmpq_poly_evaluate_fmpq(fmpq_t res, const fmpq_poly_t poly, const fmpq_t a)
{
    if (res != a)
    {
        _fmpq_poly_evaluate_fmpq(fmpq_numref(res), fmpq_denref(res),
            poly->coeffs, poly->den, poly->length, 
            fmpq_numref(a), fmpq_denref(a));
    }
    else
    {
        fmpq_t t;

        fmpq_init(t);
        fmpq_set(t, a);
        fmpq_poly_evaluate_fmpq(res, poly, t);
        fmpq_clear(t);
    }
}

