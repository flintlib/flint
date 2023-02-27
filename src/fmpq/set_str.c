/*
    Copyright (C) 2018 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int fmpq_set_str(fmpq_t x, const char *str, int base)
{
    int ans;
    mpq_t copy;

    mpq_init(copy);
    ans = mpq_set_str(copy, (char *) str, base);
    if (ans == 0)
        fmpq_set_mpq(x, copy);
    mpq_clear(copy);

    return ans;
}
