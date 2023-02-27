/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void fmpq_init_set_readonly(fmpq_t f, const mpq_t z)
{
    fmpz_init_set_readonly(fmpq_numref(f), mpq_numref(z));
    fmpz_init_set_readonly(fmpq_denref(f), mpq_denref(z));
}

