/*
    Copyright (C) 2010 Sebastian Pancratz

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

void 
_fmpq_poly_evaluate_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz * poly, 
                         const fmpz_t den, slong len, const fmpz_t a)
{
    fmpz_t d;
    
    _fmpz_poly_evaluate_horner_fmpz(rnum, poly, len, a);
    
    fmpz_init(d);
    fmpz_gcd(d, rnum, den);
    if (*d != WORD(1))
    {
        fmpz_divexact(rnum, rnum, d);
        fmpz_divexact(rden, den, d);
    }
    else
    {
        fmpz_set(rden, den);
    }
    fmpz_clear(d);
}

void 
fmpq_poly_evaluate_fmpz(fmpq_t res, const fmpq_poly_t poly, const fmpz_t a)
{
    _fmpq_poly_evaluate_fmpz(fmpq_numref(res), fmpq_denref(res), 
                             poly->coeffs, poly->den, poly->length, a);
}

