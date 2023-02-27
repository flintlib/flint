/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

double fmpq_get_d(const fmpq_t a)
{
    double d;
    mpq_t z;
    flint_mpq_init_set_readonly(z, a);
    d = mpq_get_d(z);
    flint_mpq_clear_readonly(z);
    return d;
}
