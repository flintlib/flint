/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f)
{
    flint_mpz_init_set_readonly(mpq_numref(z), fmpq_numref(f));
    flint_mpz_init_set_readonly(mpq_denref(z), fmpq_denref(f));
}

