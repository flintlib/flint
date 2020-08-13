/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void fmpq_clear_readonly(fmpq_t f)
{
    fmpz_clear_readonly(fmpq_numref(f));
    fmpz_clear_readonly(fmpq_denref(f));
}

