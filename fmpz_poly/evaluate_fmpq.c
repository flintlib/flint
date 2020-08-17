/*
    Copyright (C) 2013 Fredrik Johansson

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

void
_fmpz_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden,
                               const fmpz * f, slong len, 
                               const fmpz_t anum, const fmpz_t aden)
{
    if (len < 40 || fmpz_bits(aden) > 0.003 * len * len)
        _fmpz_poly_evaluate_horner_fmpq(rnum, rden, f, len, anum, aden);
    else
        _fmpz_poly_evaluate_divconquer_fmpq(rnum, rden, f, len, anum, aden);
}

void
fmpz_poly_evaluate_fmpq(fmpq_t res, const fmpz_poly_t f, const fmpq_t a)
{
    if (res == a)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_poly_evaluate_fmpq(t, f, a);
        fmpq_swap(res, t);
        fmpq_clear(t);
    }
    else
    {
        _fmpz_poly_evaluate_fmpq(fmpq_numref(res), fmpq_denref(res),
            f->coeffs, f->length, fmpq_numref(a), fmpq_denref(a));
    }
}

