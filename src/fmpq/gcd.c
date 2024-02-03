/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void
_fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            const fmpz_t r, const fmpz_t s)
{
   fmpz_t a, b;
   fmpz_init(a); fmpz_init(b);
   fmpz_mul(a, p, s);
   fmpz_mul(b, q, r);
   fmpz_gcd(rnum, a, b);
   fmpz_mul(rden, q, s);
   _fmpq_canonicalise(rnum, rden);
   fmpz_clear(a); fmpz_clear(b);
}

void
fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_gcd(fmpq_numref(res), fmpq_denref(res), fmpq_numref(op1),
              fmpq_denref(op1), fmpq_numref(op2), fmpq_denref(op2));
}
