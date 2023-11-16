/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "arb.h"
#include "bernoulli.h"

void
bernoulli_rev_clear(bernoulli_rev_t iter)
{
    if (iter->alloc != 0)
    {
        _fmpz_vec_clear(iter->powers, iter->alloc);
        fmpz_clear(iter->pow_error);
        arb_clear(iter->prefactor);
        arb_clear(iter->two_pi_squared);
    }
}

