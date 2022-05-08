/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"

void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
{
    fmpz_set_mpz(fmpq_numref(dest), mpq_numref(src));
    fmpz_set_mpz(fmpq_denref(dest), mpq_denref(src));
}
