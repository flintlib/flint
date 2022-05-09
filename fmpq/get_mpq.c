/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "fmpz.h"
#include "fmpq.h"

void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
{
    fmpz_get_mpz(mpq_numref(dest), fmpq_numref(src));
    fmpz_get_mpz(mpq_denref(dest), fmpq_denref(src));
}
