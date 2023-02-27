/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gmp.h"
#include "fmpq.h"

void fmpq_init_set_mpz_frac_readonly(fmpq_t z, const mpz_t num, const mpz_t den)
{
    fmpz_init_set_readonly(fmpq_numref(z), num);
    fmpz_init_set_readonly(fmpq_denref(z), den);
}
