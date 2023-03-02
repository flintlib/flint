/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

void
fmpz_lll_context_init_default(fmpz_lll_t fl)
{
    fl->delta = 0.99;
    fl->eta = 0.51;
    fl->rt = Z_BASIS;
    fl->gt = APPROX;
}
